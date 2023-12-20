#include "FLIPSimulation.h"


FLIPSimulation::FLIPSimulation()
{

}

FLIPSimulation::FLIPSimulation(float density, float width, float depth, float height, float spacing, float particleRadius, float maxParticles, float dt)
    : FluidSimulation(width, depth, height)
{
    init(density, width, depth, height, spacing, particleRadius, maxParticles, dt);
}

void FLIPSimulation::init(float density, float width, float depth, float height, float spacing, float particleRadius, float maxParticles, float dt)
{
    spacing = 1.f;
    this->dimensions = Vector3(width, depth, height);
    this->dt = dt;

    // fluid
    this->density = density;
    this->fNum = Vector3(std::floor(width / spacing) + 1, std::floor(depth / spacing) + 1, std::floor(height / spacing) + 1);

    this->cellSize = height / fNum.z; // std::max(width / this->fNum.x, height / this->fNum.y);
    this->fInvSpacing = 1.0 / spacing; //this->h;

    int nbCells = fNum.x * fNum.y * fNum.z;
    this->u = GridF(fNum);
    this->v = GridF(fNum);
    this->w = GridF(fNum);
    this->du = GridF(fNum);
    this->dv = GridF(fNum);
    this->dw = GridF(fNum);
    this->prevU = GridF(fNum);
    this->prevV = GridF(fNum);
    this->prevW = GridF(fNum);
    this->p = GridF(fNum);
    this->s = GridF(fNum);
    this->cellType = GridI(fNum);

    // particles
    this->maxParticles = maxParticles;

    this->particleDensity = GridF(fNum);
    this->particleRestDensity = 1000.0;

    this->particleRadius = particleRadius;
    this->pInvSpacing = 1.0 / (2.2 * particleRadius);
    this->pNumX = std::floor(width * this->pInvSpacing) + 1;
    this->pNumY = std::floor(depth * this->pInvSpacing) + 1;
    this->pNumZ = std::floor(height * this->pInvSpacing) + 1;
    this->pNumCells = this->pNumX * this->pNumY * this->pNumZ;

    this->numCellParticles = std::vector<int>(this->pNumCells);
    this->firstCellParticle = std::vector<int>(this->pNumCells + 1);
    this->cellParticleIds = std::vector<int>(maxParticles);

    this->particles.resize(maxParticles);
    for (int i = 0; i < maxParticles; i++) {
        Vector3 newPos(false);
        while (!newPos.isValid() || (!this->obstacleGrid.empty() && this->obstacleGrid(newPos) > 0)) {
            newPos = this->dimensions * Vector3::random(Vector3(.1f, .1f, 1.f), Vector3(.9f, .9f, 1.1f));
            newPos.z = this->dimensions.z - 3;
        }
        particles[i].position = newPos;
    }
    this->savedState = particles;
}

void FLIPSimulation::integrateParticles(double dt, const Vector3& gravity) {
    for (size_t i = 0; i < particles.size(); i++) {
        Particle& p = particles[i];
        p.velocity += gravity * dt;
        p.position = p.position + p.velocity * dt;
        p.isGhost = std::max(0, p.isGhost - 1); // New-born particle?
    }
}

void FLIPSimulation::pushParticlesApart(int numIters) {

    particleTree.build(particles);

    int numParticles = particles.size();

    for (int iter = 0; iter < numIters; iter++) {
//    #pragma omp parallel for
        for (int i = 0; i < numParticles; i++) {
            if (particles[i].isGhost > 0)
                continue;
            Vector3 p = particles[i].position;

            Vector3 displacement;
            std::vector<size_t> neighborsIDs;
            particleTree.findNeighbors(particles, particleTree.root, p, 2.f * particleRadius, neighborsIDs);
            for (auto neighborsID : neighborsIDs) {
                displacement += p - particles[neighborsID].position;
            }
            particles[i].position += displacement;
        }
    }
    return;
}

void FLIPSimulation::handleCollisions() {

    int numParticles = this->particles.size();
#pragma omp parallel for
    for (int i = 0; i < numParticles; i++) {
        Vector3& pos = particles[i].position;
        Vector3& vel = particles[i].velocity;

        Vector3 initialPoint = pos;
        Vector3 start = savedState[i].position;
        Vector3 end = pos;
        Vector3 collisionPoint(false);
        Vector3 collisionNormal(false);

        size_t last_hit_index = -1;
        int maxTries = 5;
        do {
            std::tie(collisionPoint, last_hit_index) = obstacleTriangleTree.getIntersectionAndTriangleIndex(start, end, {last_hit_index});
            if (collisionPoint.isValid()) {
                collisionNormal = obstacleTriangleTree.triangles[last_hit_index].normal;
                float distToCollision = (start - collisionPoint).norm();
                Vector3 bounce = vel.reflexion(collisionNormal).setMag((end - start).norm() - distToCollision);
                pos = collisionPoint + bounce + collisionNormal * .01f;
    //            vel *= 0.f;
                start = collisionPoint + collisionNormal * .01f;
                end = pos;
            }
            if (maxTries < 2) {
                int a = 0;
            }
            if (maxTries-- < 0) {
                pos += (pos - initialPoint) * .5f;
                break;
            }
        } while(collisionPoint.isValid());
    }
}


void FLIPSimulation::updateParticleDensity() {

    particleDensity.reset();

    int numParticles = this->particles.size();
    for (int i = 0; i < numParticles; i++) {
        Vector3 pos = particles[i].position;
        particleDensity.addValueAt(1.f, pos); // TODO : parallelisation
    }

    if (particleRestDensity == 0.0) {
        float sum = 0.0;
        int numFluidCells = 0;

        for (int i = 0; i < countCells(); i++) {
            if (cellType[i] == FLUID_CELL) {
                sum += particleDensity[i];
                numFluidCells++;
            }
        }

        if (numFluidCells > 0)
            particleRestDensity = sum / numFluidCells;
    }
    std::cout << "Total density: " << particleDensity.sum() << " (out of " << numParticles << " particles)\nMin/max density: " << particleDensity.min() << "/" << particleDensity.max() << std::endl;
}

void FLIPSimulation::transferVelocities(bool toGrid, float flipRatio)
{
    if (toGrid) {
        prevU = u;
        prevV = v;
        prevW = w;

        du.reset();
        dv.reset();
        dw.reset();
        u.reset();
        v.reset();
        w.reset();

        for (int i = 0; i < countCells(); i++) {
            cellType[i] = (s[i] == 0.0f) ? SOLID_CELL : AIR_CELL;
        }

        int numParticles = this->particles.size();
    #pragma omp parallel for
        for (int i = 0; i < numParticles; i++) {
            if (cellType(particles[i].position) == AIR_CELL)
                cellType(particles[i].position) = FLUID_CELL;
        }
    }


    int numParticles = this->particles.size();
    if (toGrid) {
        #pragma omp parallel for
        for (int i = 0; i < numParticles; i++) {
            auto& p = particles[i];

            u.addValueAt(p.velocity.x, p.position);
            v.addValueAt(p.velocity.y, p.position);
            w.addValueAt(p.velocity.z, p.position);

            du.addValueAt(1.f, p.position);
            dv.addValueAt(1.f, p.position);
            dw.addValueAt(1.f, p.position);
        }

        u.iterateParallel([&](size_t i) {
            u[i] /= (du[i] != 0 ? du[i] : 1.f);
            v[i] /= (dv[i] != 0 ? dv[i] : 1.f);
            w[i] /= (dw[i] != 0 ? dw[i] : 1.f);
        });

        // restore solid cells
        cellType.iterateParallel([&](const Vector3& pos) {
            bool solid = cellType(pos) == SOLID_CELL;
            if (solid || cellType(pos - Vector3(1, 0, 0)) == SOLID_CELL)
                u(pos) = prevU(pos);
            if (solid || cellType(pos - Vector3(0, 1, 0)) == SOLID_CELL)
                v(pos) = prevV(pos);
            if (solid || cellType(pos - Vector3(0, 0, 1)) == SOLID_CELL)
                w(pos) = prevW(pos);
        });
    } else {
#pragma omp parallel for
        for (size_t i = 0; i < numParticles; i++) {
            auto& p = particles[i];

            Vector3 vel = p.velocity;
            Vector3 picV = Vector3(u.interpolate(p.position), v.interpolate(p.position), w.interpolate(p.position)); //Vector3(getU(p.position, u), getV(p.position, v), getW(p.position, w));
            Vector3 corr = picV - Vector3(prevU.interpolate(p.position), prevV.interpolate(p.position), prevW.interpolate(p.position)); //Vector3(getU(p.position, prevU), getV(p.position, prevV), getW(p.position, prevW));
            Vector3 flipV = vel + corr;

            p.velocity = (1.f - flipRatio) * picV + flipRatio * flipV;
        }

    }
    return;
}


void FLIPSimulation::solveIncompressibility(int numIters, float dt, float overRelaxation, bool compensateDrift) {

    const float cp = density * cellSize / dt;

    for (int iter = 0; iter < numIters; iter++) {
        p.reset();
        s.reset();
        s.iterateParallel([&](int x, int y, int z) {
            if (obstacleGrid(x, y, z) > 0)
                s(x, y, z) = 0;
            else //if (x == 0 || y == 0 || z == 0 || x == s.sizeX - 1 || y == s.sizeY - 1 || z == s.sizeZ - 1)
                s(x, y, z) = 1; //s.getNumberNeighbors(x, y, z, false);
        });
        GridV3 velocities(u.getDimensions());
        velocities.iterateParallel([&](size_t i) {
            velocities[i] = Vector3(u[i], v[i], w[i]);
        });
        GridF divergences = velocities.divergence();

        cellType.iterateParallel([&](const Vector3& pos) {
            if (cellType(pos) != FLUID_CELL)
                return;

            Vector3 center = pos;
            Vector3 left = pos - Vector3(1, 0, 0);
            Vector3 right = pos + Vector3(1, 0, 0);
            Vector3 bottom = pos - Vector3(0, 1, 0);
            Vector3 top = pos + Vector3(0, 1, 0);
            Vector3 back = pos - Vector3(0, 0, 1);
            Vector3 front = pos + Vector3(0, 0, 1);

            float sumS = this->s(center) +
                         this->s(left) + this->s(right) +
                         this->s(bottom) + this->s(top) +
                         this->s(back) + this->s(front);

            if (sumS == 0.0)
                return;

            float div = divergences(pos);

            if (particleRestDensity > 0.0 && compensateDrift) {
                const float k = 1.0;
                float compression = particleDensity(center) - particleRestDensity;
                if (compression > 0.0)
                    div -= k * compression;
            }

            const float pressure = -div / sumS * overRelaxation;
            p(center) += cp * pressure;

            u(center) -= (this->s(left) * pressure);
            u(right) += (this->s(right) * pressure);
            v(center) -= (this->s(bottom) * pressure);
            v(top) += (this->s(top) * pressure);
            w(center) -= (this->s(back) * pressure);
            w(front) += (this->s(front) * pressure);
        });
//        std::cout << "Step " << iter << "\nTotal pressure: " << p.sum() << "\nMin/max pressure: " << p.min() << "/" << p.max() << std::endl;
    }
}

bool checkValidPositions(std::vector<Particle> particles) {
    bool validPositions = true;
    for (size_t i = 0; i < particles.size(); i++) {
        const auto& particle = particles[i];
        if (particle.position.z < 4.f) {
            std::cout << "Invalid #" << i << ": " << particle.position << std::endl;
            validPositions = false;
        }
    }
    if (!validPositions) {
        int a = 0;
    }
    return validPositions;
}

void FLIPSimulation::step()
{
    currentStep++;
    Vector3 gravity(0, 0, -gravityValue);
    savedState = this->particles;

    float integrationTime, pushingTime, collisionTime1, collisionTime2, densityTime, transferTime, transfer2Time, solvingTime, respawnTime, storingTime, totalTime;


    totalTime = timeIt([&]() {
        integrationTime = timeIt([&](){ integrateParticles(dt, gravity); });
        this->useVelocityForCollisionDetection = false;
        pushingTime = timeIt([&](){ pushParticlesApart(numIterations); });
        collisionTime1 = timeIt([&](){ handleCollisions(); });
        densityTime = timeIt([&](){ updateParticleDensity(); });
        transferTime = timeIt([&](){ transferVelocities(true, flipRatio); });
        solvingTime = timeIt([&](){ solveIncompressibility(numIterations, dt, overRelaxation, compensateDrift); });
        collisionTime2 = timeIt([&](){ handleCollisions(); });
        transfer2Time = timeIt([&](){ transferVelocities(false, flipRatio); });
        respawnTime = timeIt([&](){ respawnLostParticles(); });
        storingTime = timeIt([&]() { storeVelocities(); });
    });

    for (size_t i = 0; i < particles.size(); i++) {
        if (std::abs((particles[i].position - savedState[i].position).normalize().dot(Vector3(0, 0, 1))) < .9f) {
//            std::cout << "Particle #" << i << " -> " << (particles[i].position - savedState[i].position).normalize().dot(Vector3(0, 0, 1)) << std::endl;
//            particles[i] = savedState[i];
        }
    }

    /*std::cout << "Integration: " << showTime(integrationTime) << "\n"
              << "Pushing    : " << showTime(pushingTime) << "\n"
              << "Collision1 : " << showTime(collisionTime1) << "\n"
              << "Collision2 : " << showTime(collisionTime2) << "\n"
              << "Density    : " << showTime(densityTime) << "\n"
              << "Transfer   : " << showTime(transferTime) << "\n"
              << "Solving    : " << showTime(solvingTime) << "\n"
              << "Transfer 2 : " << showTime(transfer2Time) << "\n"
              << "Storage    : " << showTime(storingTime) << "\n"
              << "Total: " << showTime(totalTime) << std::endl;*/
}

void FLIPSimulation::simulate() {

    while (true) {
        step();
    }
}

void FLIPSimulation::respawnLostParticles()
{
#pragma omp parallel for
    for (size_t i = 0; i < particles.size(); i++) {
        auto& p = particles[i];
        if (!Vector3::isInBox(p.position, Vector3(), this->dimensions + Vector3(0, 0, 10000))) {
            Vector3 previousPos = p.position;
            p.position = this->dimensions * Vector3::random(Vector3(.0f, .0f, .0f), Vector3(.0f, 1.f, 1.0f));// + Vector3(0, 0, 5.1f);
            if (this->savedState.size() > i) {
                savedState[i].position = p.position - Vector3(1, 0, 0);
                p.velocity = Vector3(1, 0, 0);
                p.isGhost = 20; // countdown to the "locked" state
            }
        }
    }
}

void FLIPSimulation::storeVelocities()
{
    GridV3 currentVelocities(this->dimensions);
    GridF amount(this->dimensions);

    for (size_t i = 0; i < particles.size(); i++) {
        const auto& particle = particles[i];
        currentVelocities(particle.position) += (particle.position - savedState[i].position) / dt;
        amount(particle.position) += 1;
    }

    currentVelocities.iterateParallel([&](size_t i) {
        currentVelocities[i] = (amount[i] != 0 ? currentVelocities[i] / amount[i] : currentVelocities[i]);
    });

    if (velocitiesHistory.size() > maxAverageSize) {
        velocitiesHistory.erase(velocitiesHistory.begin(), velocitiesHistory.begin() + (velocitiesHistory.size() - maxAverageSize));
    }
    velocitiesHistory.push_back(currentVelocities);
}

GridV3 FLIPSimulation::getVelocities(int newSizeX, int newSizeY, int newSizeZ)
{
    if (_cachedStep != currentStep) {
        _cachedStep = currentStep;

        GridV3 velocities(dimensions);
        if (!velocitiesHistory.empty()) {
            for (size_t i = 0; i < std::min(averaging, int(velocitiesHistory.size())); i++) {
                const auto& velHistory = velocitiesHistory[velocitiesHistory.size() - (i + 1)];
                velocities += velHistory;
            }
            velocities /= float(velocitiesHistory.size());
        }

        this->_cachedVelocity = velocities;
    }
    return _cachedVelocity.resize(newSizeX, newSizeY, newSizeZ);
}

Vector3 FLIPSimulation::getVelocity(int x, int y, int z)
{
    return FluidSimulation::getVelocity(x, y, z);
}

void FLIPSimulation::addVelocity(int x, int y, int z, const Vector3 &amount)
{
    float radiusEffect = 1.f;
    Vector3 effectArea(x, y, z);
    for (auto& particle : particles) {
        if ((particle.position - effectArea).norm2() < radiusEffect * radiusEffect)
            particle.velocity += amount;
    }
}

void FLIPSimulation::reset()
{
    this->init(this->density, this->dimensions.x, this->dimensions.y, this->dimensions.z, 1.f / this->fInvSpacing, this->particleRadius, this->maxParticles, this->dt);
}

float FLIPSimulation::getU(const Vector3 &p, const GridF &uGrid) const
{
    Vector3 pp = toUcoords(p);
    return uGrid.interpolate(pp);
    /*
    Vector3 pp = p;
    pp.x -= 1.5f;
    Vector3 pp2 = p;
    pp2.x += 1.5f;
    float t = pp.x - int(pp.x);
    return uGrid(pp) * (1.f - t) + uGrid(pp2) * t;*/
}

float FLIPSimulation::getV(const Vector3 &p, const GridF &vGrid) const
{
    Vector3 pp = toVcoords(p);
    return vGrid.interpolate(pp);
    /*
    Vector3 pp = p;
    pp.y -= 1.5f;
    Vector3 pp2 = p;
    pp2.y += 1.5f;
    float t = pp.y - int(pp.y);
    return vGrid(pp) * (1.f - t) + vGrid(pp2) * t;
    */
}

float FLIPSimulation::getW(const Vector3 &p, const GridF &wGrid) const
{
    Vector3 pp = toWcoords(p);
    return wGrid.interpolate(pp);
    /*
    Vector3 pp = p;
    pp.z -= 1.5f;
    Vector3 pp2 = p;
    pp2.z += 1.5f;
    float t = pp.z - int(pp.z);
    return wGrid(pp) * (1.f - t) + wGrid(pp2) * t;
    */
}
