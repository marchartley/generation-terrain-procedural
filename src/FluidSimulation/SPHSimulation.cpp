#include "SPHSimulation.h"

using namespace SPH;

SPHSimulation::SPHSimulation()
{
    this->dimensions = Vector3(100, 100, 30);
    this->nbParticles = 50000;
    this->dt = 0.01f;
    this->damping = .99f;
}

SPHSimulation::~SPHSimulation()
{
    delete this->tree;
}

void SPHSimulation::computeNeighbors() {
    size_t nbParticles = particles.size();
    precomputedNeighbors.resize(nbParticles);
#pragma omp parallel for
    for (size_t i = 0; i < nbParticles; i++) {
        Particle& particle = particles[i];
        precomputedNeighbors[i] = getNeighbors(particle.position, particle.smoothingRadius);
        if (!Vector3::isInBox(particle.position, Vector3(), this->dimensions)) {
            particle.isGhost = true;
            particle.velocity *= 0.f;
            particle.force *= 0.f;
        }
    }
}
void SPHSimulation::initialize(std::vector<std::vector<Vector3> > meshBoundaries) {
    t = 0.f;
    computeTime = 0.f;
    this->particles.resize(nbParticles);

    Particle particleModel;
    particleModel.velocity = Vector3(0.0f, 0.0f, 0.0f);
    particleModel.density = 0.0f;
    particleModel.pressure = 0.0f;
    particleModel.mass = 10.f;
    particleModel.smoothingRadius = 1.f;
    particleModel.gasConstant = 5.f;
    particleModel.restDensity = 997.f;
    particleModel.viscosity = .0001f;
    particleModel.isGhost = false;

    for (int i = 0; i < nbParticles; i++) {
        Particle particle = particleModel;
        particle.position = Vector3::random(Vector3(1, 1, 1), dimensions - Vector3(1, 1, 1));
        particles[i] = particle;
//        cells.at(particle.position).push_back(i);
    }
/*
    // Add ghost particles at the boundaries
    Vector3 nbGhosts = dimensions / (particleModel.smoothingRadius * .5f);
    int nbX = nbGhosts.x;
    int nbY = nbGhosts.y;
    int nbZ = nbGhosts.z;
    for (int x = 0; x <= nbX; x++) {
        for (int y = 0; y <= nbY; y++) {
            for (int z = 0; z <= nbZ; z++) {
                if (x == 0 || y == 0 || z == 0 || x == nbX || y == nbY || z == nbZ) {
                    Particle ghostParticle = particleModel;
                    ghostParticle.position = Vector3(x, y, z) * .5f;
                    ghostParticle.isGhost = true;

                    particles.push_back(ghostParticle);
//                    cells.at(ghostParticle.position).push_back(particles.size() - 1);
                }
            }
        }
    }*/

    // Add ghost particles on the mesh boundaries
    for (const auto& triangle : meshBoundaries) {
        Vector3 v1 = triangle[0];
        Vector3 v2 = triangle[1];
        Vector3 v3 = triangle[2];

        if (!Vector3::isInBox(v1, Vector3(), this->dimensions) && !Vector3::isInBox(v2, Vector3(), this->dimensions) && !Vector3::isInBox(v3, Vector3(), this->dimensions))
            continue;

        // Calculate the edge vectors
        Vector3 edge1 = v2 - v1;
        Vector3 edge2 = v3 - v1;

        // Calculate the number of particles to place along each edge
        int numParticlesEdge1 = std::ceil(edge1.magnitude() / particleModel.smoothingRadius);
        int numParticlesEdge2 = std::ceil(edge2.magnitude() / particleModel.smoothingRadius);

        // Place particles along the edges and the face of the triangle
        for (int i = 0; i <= numParticlesEdge1; ++i) {
            for (int j = 0; j <= numParticlesEdge2; ++j) {
                if (i + j <= std::max(numParticlesEdge1, numParticlesEdge2)) {
                    Particle ghostParticle = particleModel;
                    ghostParticle.position = v1 + edge1 *((i + 1.f) / (numParticlesEdge1 + 1.f)) + edge2 * ((j + 1.f) / (numParticlesEdge2 + 1.f));
                    ghostParticle.isGhost = true;
                    particles.push_back(ghostParticle);
//                    cells.at(ghostParticle.position).push_back(particles.size() - 1);
                }
            }
        }
    }
    for (int i = 0; i < int(particles.size()); i++)
        particles[i].index = i;
    tree = new KDTree(particles);

    std::cout << nbParticles << " fluid particles for " << particles.size() - nbParticles << " ghost particles" << std::endl;
}

void SPHSimulation::computeDensityAndPressure() {
    float B = 2.2e9;  // Bulk modulus
    float gamma = 7.0;  // Polytropic exponent

    size_t nbParticles = particles.size();
#pragma omp parallel for
    for (size_t i = 0; i < nbParticles; i++) {
        Particle& particle = particles[i];
//        if (particle.isGhost) continue;
        particle.density = 0.0f;
        for (size_t j : precomputedNeighbors[i]) { //getNeighbors(particle.position, particle.smoothingRadius)) {
            if (i == j)
                continue;
            Particle& neighbor = particles[j];
            Vector3 diff = particle.position - neighbor.position;
            float r2 = diff.dot(diff);
            if (r2 < particle.smoothingRadius * particle.smoothingRadius) {
                float r = std::sqrt(r2);
                float W = particle.mass * (315.0f / (64.0f * M_PI * std::pow(particle.smoothingRadius, 9))) * std::pow(particle.smoothingRadius * particle.smoothingRadius - r * r, 3);
                particle.density += W;
            }
        }
        particle.pressure = B * (std::pow(particle.density / particle.restDensity, gamma) - 1);

        particle.pressure = std::max(0.1f, particle.pressure);
        particle.density = std::max(0.1f, particle.density);
    }
}

void SPHSimulation::computeForces() {
    Vector3 gravity = Vector3(0.0f, 0.0f, -9.8f);
    float artificialPressure = 0.1f;  // Artificial pressure coefficient

    size_t nbParticles = particles.size();
#pragma omp parallel for
    for (size_t i = 0; i < nbParticles; i++) {
        Particle& particle = particles[i];
        if (particle.isGhost) continue;
        Vector3 pressureForce(0.0f, 0.0f, 0.0f);
        Vector3 viscosityForce(0.0f, 0.0f, 0.0f);

        for (size_t j : precomputedNeighbors[i]) { //getNeighbors(particle.position, particle.smoothingRadius)) {
            if (i == j)
                continue;
            Particle& neighbor = particles[j];
            Vector3 diff = particle.position - neighbor.position;
            float r2 = diff.dot(diff);
            if (r2 == 0) {
                pressureForce += Vector3::random();
            } else if (r2 < particle.smoothingRadius * particle.smoothingRadius) {
                float r = std::sqrt(r2);
                float Wp = particle.mass * (-45.0f / (M_PI * std::pow(particle.smoothingRadius, 6))) * std::pow(particle.smoothingRadius - r, 2);
                pressureForce += (particle.pressure + neighbor.pressure) / (2.0f * neighbor.density) * Wp * diff / r;
                float Wv = particle.mass * (45.0f / (M_PI * std::pow(particle.smoothingRadius, 6))) * (particle.smoothingRadius - r);
                viscosityForce += particle.viscosity * (neighbor.velocity - particle.velocity) / neighbor.density * Wv;
            }
        }
        // Add artificial pressure term
        pressureForce += artificialPressure * particle.density * particle.density * particle.pressure;
        particle.force = pressureForce + viscosityForce + gravity * particle.density;
    }
}

void SPHSimulation::integrate() {
    // Adaptative time step using CLF conditions
    float maxSpeed = 0.f;
    for (Particle& particle : particles) {
        maxSpeed = std::max(maxSpeed, particle.velocity.norm2());
    }

    dt = 0.1f;
    if (maxSpeed > 0.f) {
        dt = clamp(0.1f * particles[0].smoothingRadius / std::sqrt(maxSpeed), 0.001f, dt);
    }

    size_t nbParticles = particles.size();
#pragma omp parallel for
    for (size_t i = 0; i < nbParticles; i++) {
        Particle& particle = particles[i];
        if (particle.isGhost) continue;
//        particle.force += Vector3(1, 0, 0);
        Vector3 acceleration = particle.force / particle.density;
        Vector3 velocityHalfStep = (particle.velocity + 0.5f * dt * acceleration);
        particle.position += dt * velocityHalfStep;
        Vector3 newAcceleration = particle.force / particle.density;
        particle.velocity = (damping * (velocityHalfStep + 0.5f * dt * newAcceleration));
    }
    this->t += dt;
}

void SPHSimulation::handleCollisions() {
    float restitution = .99f; // 0.5f;

    size_t nbParticles = particles.size();
#pragma omp parallel for
    for (size_t i = 0; i < nbParticles; i++) {
        Particle& particle = particles[i];
        if (particle.isGhost) continue;
        Vector3 displacement = Vector3(0.0f, 0.0f, 0.0f);
        // Check for collisions with the ghost particles
        for (size_t j : precomputedNeighbors[i]) { //getNeighbors(particle.position, particle.smoothingRadius)) {
            Particle& neighbor = particles[j];
            if (neighbor.isGhost) {
                Vector3 diff = particle.position - neighbor.position;
                float r2 = diff.dot(diff);
                if (r2 < particle.smoothingRadius * particle.smoothingRadius) {
                    float r = std::sqrt(r2);
                    displacement += diff.normalized() * (particle.smoothingRadius - r); //(particle.smoothingRadius - r) * diff / r;
                    particle.velocity *= -(restitution); // * (r / particle.smoothingRadius));
                }
            }
        }
        // Move the particle away from the ghost particles
        particle.position += displacement;
    }
}

std::vector<size_t> SPHSimulation::getNeighbors(Vector3& position, float maxDistance) {
    std::vector<size_t> neighborsIndices;
    tree->findNeighbors(particles, tree->root, position, maxDistance * 2.f, neighborsIndices);
    return neighborsIndices;
}
void SPHSimulation::relaxDensity() {
    float threshold = 0.1f;  // Set this to a suitable value
    int maxIterations = 2;  // Set this to a suitable value

    for (int iteration = 0; iteration < maxIterations; ++iteration) {
        // Compute the density error for each particle
        std::vector<float> densityErrors(particles.size());
        float maxDensityError = 0.0f;
        for (size_t i = 0; i < particles.size(); ++i) {
            Particle& particle = particles[i];
            float densityError = particle.density - particle.restDensity;
            densityErrors[i] = densityError;
            maxDensityError = std::max(maxDensityError, std::abs(densityError));
        }

        // If the maximum density error is below the threshold, we're done
        if (maxDensityError < threshold) {
            break;
        }

        // Otherwise, adjust the pressures of the particles to correct the density error
        for (size_t i = 0; i < particles.size(); ++i) {
            Particle& particle = particles[i];
            particle.pressure += particle.gasConstant * densityErrors[i];
        }

        // Recompute the forces and integrate
        computeForces();
        integrate();
    }
}


void SPHSimulation::step() {
    float timeKDTree = timeIt([=]() { this->computeNeighbors(); });
    std::vector<Particle> previousState = particles;
    bool nans = false;
    for (const auto& p : particles) nans |= (p.position.x != p.position.x);
    float timeDensity = timeIt([=]() { this->computeDensityAndPressure(); });
    for (const auto& p : particles) nans |= (p.position.x != p.position.x);
    float timeForces = timeIt([=]() { this->computeForces(); });
    for (const auto& p : particles) nans |= (p.position.x != p.position.x);
    float timeRelaxation = timeIt([=]() { this->relaxDensity(); });
    for (const auto& p : particles) nans |= (p.position.x != p.position.x);
    float timeIntegration = timeIt([=]() { this->integrate(); });
    for (const auto& p : particles) nans |= (p.position.x != p.position.x);
    float timeCollisions = timeIt([=]() { this->handleCollisions(); });
    for (const auto& p : particles) nans |= (p.position.x != p.position.x);
    computeTime += (timeKDTree + timeDensity + timeRelaxation + timeForces + timeIntegration + timeCollisions);
    std::cout << "KD-Tree: " << timeKDTree << "ms - Density: " << timeDensity << "ms - Relaxation: " << timeRelaxation << "ms - Forces: " << timeForces << "ms - Integration: " << timeIntegration << "ms - Collisions: " << timeCollisions << "ms ... t = " << t << " (" << computeTime/1000.f << "s)" << std::endl;

    if (nans) {
        std::cout << "NaN found!" << std::endl;
        particles = previousState;
    }
}

Matrix3<Vector3> SPHSimulation::getVelocities(int newSizeX, int newSizeY, int newSizeZ)
{
    Matrix3<Vector3> velocities(sizeX, sizeY, sizeZ);
    Matrix3<float> amount(sizeX, sizeY, sizeZ);

    for (auto& particle : this->particles) {
        velocities[particle.position] += particle.velocity;
        amount[particle.position] += 1;
    }
    for (size_t i = 0; i < velocities.size(); i++)
        velocities[i] /= amount[i];
    return velocities.resize(newSizeX, newSizeY, newSizeZ);
}

void SPHSimulation::addVelocity(int x, int y, int z, Vector3 amount)
{
    float radiusEffect = 1.f;
    Vector3 effectArea(x, y, z);
    for (auto& particleIdx : this->getNeighbors(effectArea, radiusEffect)) {
        this->particles[particleIdx].velocity += amount;
    }
}

KDNode::KDNode(size_t pIndex, int a) : particleIndex(pIndex), left(NULL), right(NULL), axis(a) {}


KDTree::KDTree(std::vector<Particle> &particles) {
    root = build(particles, 0);
}

KDNode* KDTree::build(std::vector<Particle> particles, int depth) {
    if (particles.empty()) {
        return NULL;
    }

    int axis = depth % 3;
    std::sort(particles.begin(), particles.end(), [axis](Particle a, Particle b) {
        return a.position[axis] < b.position[axis];
    });

    int median = particles.size() / 2;
    KDNode* node = new KDNode(particles[median].index, axis);

    std::vector<Particle> left(particles.begin(), particles.begin() + median);
    std::vector<Particle> right(particles.begin() + median + 1, particles.end());

    node->left = build(left, depth + 1);
    node->right = build(right, depth + 1);

    return node;
}

void KDTree::findNeighbors(std::vector<Particle> &particles, KDNode *node, Vector3 &position, float maxDistance, std::vector<size_t> &neighbors) {
    if (node == NULL) {
        return;
    }

    float d = position[node->axis] - particles[node->particleIndex].position[node->axis];
    KDNode* nearest = d < 0 ? node->left : node->right;
    KDNode* farthest = d < 0 ? node->right : node->left;

    findNeighbors(particles, nearest, position, maxDistance, neighbors);

    if (d * d < maxDistance * maxDistance) {
        if ((particles[node->particleIndex].position - position).magnitude() < maxDistance) {
            neighbors.push_back(node->particleIndex);
        }

        findNeighbors(particles, farthest, position, maxDistance, neighbors);
    }
}
