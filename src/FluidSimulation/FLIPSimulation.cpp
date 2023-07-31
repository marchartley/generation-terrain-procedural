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
    this->fNumX = std::floor(width / spacing) + 1;
    this->fNumY = std::floor(depth / spacing) + 1;
    this->fNumZ = std::floor(height / spacing) + 1;
    this->h = height / fNumZ; // std::max(width / this->fNumX, height / this->fNumY);
    this->fInvSpacing = 1.0 / this->h;
    this->fNumCells = this->fNumX * this->fNumY * this->fNumZ;

    this->u = std::vector<float>(this->fNumCells);
    this->v = std::vector<float>(this->fNumCells);
    this->w = std::vector<float>(this->fNumCells);
    this->du = std::vector<float>(this->fNumCells);
    this->dv = std::vector<float>(this->fNumCells);
    this->dw = std::vector<float>(this->fNumCells);
    this->prevU = std::vector<float>(this->fNumCells);
    this->prevV = std::vector<float>(this->fNumCells);
    this->prevW = std::vector<float>(this->fNumCells);
    this->p = std::vector<float>(this->fNumCells);
    this->s = std::vector<float>(this->fNumCells);
    this->cellType = std::vector<int>(this->fNumCells, FLUID_CELL);
    this->cellColor = std::vector<float>(3 * this->fNumCells);

    // particles
    this->maxParticles = maxParticles;

    this->particleDensity = std::vector<float>(this->fNumCells);
    this->particleRestDensity = 0.0;

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
        particles[i].position = Vector3::random(this->dimensions * Vector3(.5f, 1.f, 1.f)) + Vector3(0, 0, this->dimensions.z);
    }

    //    this->numParticles = 0;
}

void FLIPSimulation::integrateParticles(double dt, Vector3 gravity) {
    for (Particle& p : particles) {
        p.velocity += gravity * dt;
        p.position = Vector3::wrap(p.position + p.velocity * dt, Vector3(), this->dimensions);
    }
}

void FLIPSimulation::pushParticlesApart(int numIters) {
//    float colorDiffusionCoeff = 0.001;

    // Count particles per cell
    std::fill(numCellParticles.begin(), numCellParticles.end(), 0);

    int numParticles = this->particles.size();
#pragma omp parallel for
    for (int i = 0; i < numParticles; i++) {
        Vector3 pos = particles[i].position;

        int xi = clamp((int)std::floor(pos.x * pInvSpacing), 0, pNumX - 1);
        int yi = clamp((int)std::floor(pos.y * pInvSpacing), 0, pNumY - 1);
        int zi = clamp((int)std::floor(pos.z * pInvSpacing), 0, pNumZ - 1);
        int cellNr = (xi * pNumY + yi) * pNumZ + zi;
        numCellParticles[cellNr]++;
    }

    // Partial sums
    int first = 0;
    for (int i = 0; i < pNumCells; i++) {
        first += numCellParticles[i];
        firstCellParticle[i] = first;
    }
    firstCellParticle[pNumCells] = first; // guard

    // Fill particles into cells
#pragma omp parallel for
    for (int i = 0; i < numParticles; i++) {
        Vector3 pos = particles[i].position;

        int xi = clamp((int)std::floor(pos.x * pInvSpacing), 0, pNumX - 1);
        int yi = clamp((int)std::floor(pos.y * pInvSpacing), 0, pNumY - 1);
        int zi = clamp((int)std::floor(pos.z * pInvSpacing), 0, pNumZ - 1);
        int cellNr = (xi * pNumY + yi) * pNumZ + zi;

        firstCellParticle[cellNr]--;
        cellParticleIds[firstCellParticle[cellNr]] = i;
    }

    // Push particles apart
    float minDist = 2.0 * particleRadius;
    float minDist2 = minDist * minDist;

    std::vector<bool> lockedParticles(numParticles);
    Vector3 radiusBox = Vector3(1, 1, 1) * particleRadius;
//#pragma omp parallel for
//    for (int i = 0; i < numParticles; i++) {
//        auto closeBorders = obstacleTrianglesOctree->queryRange(particles[i].position - radiusBox, particles[i].position + radiusBox);
//        lockedParticles[i] = closeBorders.size() > 0;
//    }

    for (int iter = 0; iter < numIters; iter++) {
    #pragma omp parallel for
        for (int i = 0; i < numParticles; i++) {
            Vector3 p = particles[i].position;

            int pxi = std::floor(p.x * pInvSpacing);
            int pyi = std::floor(p.y * pInvSpacing);
            int pzi = std::floor(p.z * pInvSpacing);

            int x0 = std::max(pxi - 1, 0);
            int y0 = std::max(pyi - 1, 0);
            int z0 = std::max(pzi - 1, 0);
            int x1 = std::min(pxi + 1, pNumX - 1);
            int y1 = std::min(pyi + 1, pNumY - 1);
            int z1 = std::min(pzi + 1, pNumZ - 1);

            for (int xi = x0; xi <= x1; xi++) {
                for (int yi = y0; yi <= y1; yi++) {
                    for (int zi = z0; zi <= z1; zi++) {
                        int cellNr = (xi * pNumY + yi) * pNumZ + zi;
                        int first = firstCellParticle[cellNr];
                        int last = firstCellParticle[cellNr + 1];

                        for (int j = first; j < last; j++) {
                            int id = cellParticleIds[j];
                            if (id == i)
                                continue;

                            Vector3 q = particles[id].position;
                            Vector3 d = q - p;
                            float d2 = d.dot(d);

                            if (d2 > minDist2 || d2 == 0.0)
                                continue;

                            float dLen = std::sqrt(d2);
                            float s = 0.5 * (minDist - dLen) / dLen;
                            d *= s;
//                            if (!lockedParticles[i])
                                particles[i].position -= d;
//                            if (!lockedParticles[id])
                                particles[id].position += d;

                            // Diffuse colors
//                            for (int k = 0; k < 3; k++) {
//                                float color0 = particleColor[i][k];
//                                float color1 = particleColor[id][k];
//                                float color = (color0 + color1) * 0.5;
//                                particleColor[i][k] = color0 + (color - color0) * colorDiffusionCoeff;
//                                particleColor[id][k] = color1 + (color - color1) * colorDiffusionCoeff;
//                            }
                        }
                    }
                }
            }
        }
    }
}

void FLIPSimulation::handleCollisions() {
//    float h = 1.0 / fInvSpacing;
//    float r = particleRadius;
//    float or2 = obstacleRadius * obstacleRadius;
//    float minDist = obstacleRadius + r;
//    float minDist2 = minDist * minDist;

//    Vector3 min(h + r, h + r, h + r);
//    Vector3 max((fNumX - 1) * h - r, (fNumY - 1) * h - r, (fNumZ - 1) * h - r);

    int numParticles = this->particles.size();

#pragma omp parallel for
    for (int i = 0; i < numParticles; i++) {
        Vector3& pos = particles[i].position;
        Vector3& vel = particles[i].velocity;

        Vector3 startPos = (useVelocityForCollisionDetection ? pos : savedState[i].position);
        Vector3 endPos = (useVelocityForCollisionDetection ? pos + vel.normalized() : pos);
//        Vector3 diff = endPos - startPos;
        auto [collisionPoint, collisionNormal] = obstacleTriangleTree.getIntersectionAndNormal(startPos, endPos);
        if (collisionPoint.isValid()) {
            if (useVelocityForCollisionDetection) {
                vel = vel.normalized().reflexion(collisionNormal) * vel.norm(); //collisionPoint - startPos;
                if (vel.dot(collisionNormal) < 0)
                    vel *= -1.f;
                pos = collisionPoint + vel;
            } else {
                pos = collisionPoint + (startPos - collisionPoint) * .1f;
            }
        }
        /*std::vector<OctreeNodeData> nearbyTriangles = obstacleTrianglesOctree->queryRange(startPos - diff * 3.f, endPos + diff * 3.f);
        // Check for intersections with nearby triangles
        for (auto& triangleData : nearbyTriangles) {
            auto& triangle = this->triangles[triangleData.index];
            Vector3 collisionPoint = Collision::segmentToTriangleCollision(startPos - diff, endPos + diff, triangle[0], triangle[1], triangle[2]);
            if (collisionPoint.isValid()) {
                if (useVelocityForCollisionDetection) {
                    Vector3 normal = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]).normalize();
                    vel = vel.normalized().reflexion(normal) * vel.norm(); //collisionPoint - startPos;
                    if (vel.dot(normal) < 0)
                        vel *= -1.f;
                    pos = collisionPoint + vel;
                } else {
                    pos = collisionPoint + (startPos - collisionPoint) * .1f;
                }
                break;
            }
        }*/
    }
}

void FLIPSimulation::updateParticleDensity() {
    int n = fNumY;
    int m = fNumZ;
    float h = this->h;
    float h1 = fInvSpacing;
    float h2 = 0.5 * h;

    std::fill(particleDensity.begin(), particleDensity.end(), 0.0);

    int numParticles = this->particles.size();
#pragma omp parallel for
    for (int i = 0; i < numParticles; i++) {
        Vector3 pos = particles[i].position;

        pos.x = clamp(pos.x, h, (fNumX - 1) * h);
        pos.y = clamp(pos.y, h, (fNumY - 1) * h);
        pos.z = clamp(pos.z, h, (fNumZ - 1) * h);

        int x0 = floor((pos.x - h2) * h1);
        float tx = ((pos.x - h2) - x0 * h) * h1;
        int x1 = std::min(x0 + 1, fNumX-2);

        int y0 = floor((pos.y - h2) * h1);
        float ty = ((pos.y - h2) - y0 * h) * h1;
        int y1 = std::min(y0 + 1, fNumY-2);

        int z0 = floor((pos.z - h2) * h1);
        float tz = ((pos.z - h2) - z0 * h) * h1;
        int z1 = std::min(z0 + 1, fNumZ-2);

        float sx = 1.0 - tx;
        float sy = 1.0 - ty;
        float sz = 1.0 - tz;

        if (x0 < fNumX && y0 < fNumY && z0 < fNumZ)
            particleDensity[x0 * n * m + y0 * m + z0] += sx * sy * sz;
        if (x1 < fNumX && y0 < fNumY && z0 < fNumZ)
            particleDensity[x1 * n * m + y0 * m + z0] += tx * sy * sz;
        if (x1 < fNumX && y1 < fNumY && z0 < fNumZ)
            particleDensity[x1 * n * m + y1 * m + z0] += tx * ty * sz;
        if (x0 < fNumX && y1 < fNumY && z0 < fNumZ)
            particleDensity[x0 * n * m + y1 * m + z0] += sx * ty * sz;
        if (x0 < fNumX && y0 < fNumY && z1 < fNumZ)
            particleDensity[x0 * n * m + y0 * m + z1] += sx * sy * tz;
        if (x1 < fNumX && y0 < fNumY && z1 < fNumZ)
            particleDensity[x1 * n * m + y0 * m + z1] += tx * sy * tz;
        if (x1 < fNumX && y1 < fNumY && z1 < fNumZ)
            particleDensity[x1 * n * m + y1 * m + z1] += tx * ty * tz;
        if (x0 < fNumX && y1 < fNumY && z1 < fNumZ)
            particleDensity[x0 * n * m + y1 * m + z1] += sx * ty * tz;
    }

    if (particleRestDensity == 0.0) {
        float sum = 0.0;
        int numFluidCells = 0;

        for (int i = 0; i < fNumCells; i++) {
            if (cellType[i] == FLUID_CELL) {
                sum += particleDensity[i];
                numFluidCells++;
            }
        }

        if (numFluidCells > 0)
            particleRestDensity = sum / numFluidCells;
    }
}
void FLIPSimulation::transferVelocities(bool toGrid, float flipRatio)
{
    int n = fNumY;
    float h = this->h;
    float h1 = this->fInvSpacing;
    float h2 = 0.5f * h;

    if (toGrid) {
        prevU = u;
        prevV = v;

        du.assign(fNumCells, 0.0f);
        dv.assign(fNumCells, 0.0f);
        u.assign(fNumCells, 0.0f);
        v.assign(fNumCells, 0.0f);

        for (int i = 0; i < fNumCells; i++) {
            cellType[i] = (s[i] == 0.0f) ? SOLID_CELL : AIR_CELL;
        }

        int numParticles = this->particles.size();
    #pragma omp parallel for
        for (int i = 0; i < numParticles; i++) {
            float x = particles[i].position.x;
            float y = particles[i].position.y;
            float z = particles[i].position.z;
            int xi = clamp((int)(x * h1), 0, fNumX - 1);
            int yi = clamp((int)(y * h1), 0, fNumY - 1);
            int zi = clamp((int)(z * h1), 0, fNumZ - 1);
            int cellNr = xi * n * fNumZ + yi * fNumZ + zi;
            if (cellType[cellNr] == AIR_CELL) {
                cellType[cellNr] = FLUID_CELL;
            }
        }
    }

    for (int component = 0; component < 3; component++) {
        float dx = (component == 0) ? 0.0f : h2;
        float dy = (component == 1) ? 0.0f : h2;
        float dz = (component == 2) ? 0.0f : h2;

        std::vector<float>& f = (component == 0) ? u : ((component == 1) ? v : w);
        std::vector<float>& prevF = (component == 0) ? prevU : ((component == 1) ? prevV : prevW);
        std::vector<float>& d = (component == 0) ? du : ((component == 1) ? dv : dw);

        int numParticles = this->particles.size();
#pragma omp parallel for
        for (int i = 0; i < numParticles; i++) {
            float x = particles[i].position.x;
            float y = particles[i].position.y;
            float z = particles[i].position.z;

            x = clamp(x, h, (fNumX - 1) * h);
            y = clamp(y, h, (fNumY - 1) * h);
            z = clamp(z, h, (fNumZ - 1) * h);

            int x0 = std::min((int)((x - dx) * h1), fNumX - 2);
            float tx = (x - dx - x0 * h) * h1;
            int x1 = std::min(x0 + 1, fNumX - 2);

            int y0 = std::min((int)((y - dy) * h1), fNumY - 2);
            float ty = (y - dy - y0 * h) * h1;
            int y1 = std::min(y0 + 1, fNumY - 2);

            int z0 = std::min((int)((z - dz) * h1), fNumZ - 2);
            float tz = (z - dz - z0 * h) * h1;
            int z1 = std::min(z0 + 1, fNumZ - 2);

            float sx = 1.0f - tx;
            float sy = 1.0f - ty;
            float sz = 1.0f - tz;

            float d0 = sx * sy * sz;
            float d1 = tx * sy * sz;
            float d2 = tx * ty * sz;
            float d3 = sx * ty * sz;
            float d4 = sx * sy * tz;
            float d5 = tx * sy * tz;
            float d6 = tx * ty * tz;
            float d7 = sx * ty * tz;

            int nr0 = x0 * n * fNumZ + y0 * fNumZ + z0;
            int nr1 = x1 * n * fNumZ + y0 * fNumZ + z0;
            int nr2 = x1 * n * fNumZ + y1 * fNumZ + z0;
            int nr3 = x0 * n * fNumZ + y1 * fNumZ + z0;
            int nr4 = x0 * n * fNumZ + y0 * fNumZ + z1;
            int nr5 = x1 * n * fNumZ + y0 * fNumZ + z1;
            int nr6 = x1 * n * fNumZ + y1 * fNumZ + z1;
            int nr7 = x0 * n * fNumZ + y1 * fNumZ + z1;

            if (toGrid) {
                float pv = particles[i].velocity[component];
                f[nr0] += pv * d0;
                f[nr1] += pv * d1;
                f[nr2] += pv * d2;
                f[nr3] += pv * d3;
                f[nr4] += pv * d4;
                f[nr5] += pv * d5;
                f[nr6] += pv * d6;
                f[nr7] += pv * d7;

                d[nr0] += d0;
                d[nr1] += d1;
                d[nr2] += d2;
                d[nr3] += d3;
                d[nr4] += d4;
                d[nr5] += d5;
                d[nr6] += d6;
                d[nr7] += d7;
            }
            else {
                int offset = (component == 0) ? n * fNumZ : (component == 1) ? fNumZ : 1;
                float valid0 = (cellType[nr0] != AIR_CELL || cellType[nr0 - offset] != AIR_CELL) ? 1.0f : 0.0f;
                float valid1 = (cellType[nr1] != AIR_CELL || cellType[nr1 - offset] != AIR_CELL) ? 1.0f : 0.0f;
                float valid2 = (cellType[nr2] != AIR_CELL || cellType[nr2 - offset] != AIR_CELL) ? 1.0f : 0.0f;
                float valid3 = (cellType[nr3] != AIR_CELL || cellType[nr3 - offset] != AIR_CELL) ? 1.0f : 0.0f;
                float valid4 = (cellType[nr4] != AIR_CELL || cellType[nr4 - offset] != AIR_CELL) ? 1.0f : 0.0f;
                float valid5 = (cellType[nr5] != AIR_CELL || cellType[nr5 - offset] != AIR_CELL) ? 1.0f : 0.0f;
                float valid6 = (cellType[nr6] != AIR_CELL || cellType[nr6 - offset] != AIR_CELL) ? 1.0f : 0.0f;
                float valid7 = (cellType[nr7] != AIR_CELL || cellType[nr7 - offset] != AIR_CELL) ? 1.0f : 0.0f;

                float v = particles[i].velocity[component];
                float totalValid = valid0 + valid1 + valid2 + valid3 + valid4 + valid5 + valid6 + valid7;
                float invTotalValid = 1.0f / totalValid;

                if (totalValid > 0.0f) {
                    float picV = (valid0 * d0 * f[nr0] + valid1 * d1 * f[nr1] + valid2 * d2 * f[nr2] + valid3 * d3 * f[nr3] +
                                  valid4 * d4 * f[nr4] + valid5 * d5 * f[nr5] + valid6 * d6 * f[nr6] + valid7 * d7 * f[nr7]) * invTotalValid;
                    float corr = (valid0 * d0 * (f[nr0] - prevF[nr0]) + valid1 * d1 * (f[nr1] - prevF[nr1]) +
                                  valid2 * d2 * (f[nr2] - prevF[nr2]) + valid3 * d3 * (f[nr3] - prevF[nr3]) +
                                  valid4 * d4 * (f[nr4] - prevF[nr4]) + valid5 * d5 * (f[nr5] - prevF[nr5]) +
                                  valid6 * d6 * (f[nr6] - prevF[nr6]) + valid7 * d7 * (f[nr7] - prevF[nr7])) * invTotalValid;
                    float flipV = v + corr;

                    particles[i].velocity[component] = (1.0f - flipRatio) * picV + flipRatio * flipV;
                }
            }
        }

        if (toGrid) {
            for (int i = 0; i < fNumCells; i++) {
                if (d[i] > 0.0f) {
                    f[i] /= d[i];
                }
            }

            // restore solid cells
            #pragma omp parallel for collapse(3)
            for (int i = 0; i < fNumX; i++) {
                for (int j = 0; j < fNumY; j++) {
                    for (int k = 0; k < fNumZ; k++) {
                        int idx = i * n * fNumZ + j * fNumZ + k;
                        bool solid = (cellType[idx] == SOLID_CELL);
                        if (solid || (i > 0 && cellType[(i - 1) * n * fNumZ + j * fNumZ + k] == SOLID_CELL)) {
                            u[idx] = prevU[idx];
                        }
                        if (solid || (j > 0 && cellType[i * n * fNumZ + (j - 1) * fNumZ + k] == SOLID_CELL)) {
                            v[idx] = prevV[idx];
                        }
                        if (solid || (k > 0 && cellType[i * n * fNumZ + j * fNumZ + (k - 1)] == SOLID_CELL)) {
                            w[idx] = prevW[idx];
                        }
                    }
                }
            }
        }
    }
}


void FLIPSimulation::solveIncompressibility(int numIters, float dt, float overRelaxation, bool compensateDrift) {
    std::fill(p.begin(), p.end(), 0.0);
    prevU = u;
    prevV = v;
    prevW = w;

    int n = fNumY;
    int m = fNumZ;
    float cp = density * h / dt;

    // The loop over fNumCells is redundant here as u, v, w are not being updated.

    for (int iter = 0; iter < numIters; iter++) {

        for (int i = 1; i < fNumX - 1; i++) {
            for (int j = 1; j < fNumY - 1; j++) {
                for (int k = 1; k < fNumZ - 1; k++) {

                    if (cellType[(i * n + j) * m + k] != FLUID_CELL)
                        continue;

                    int center = (i * n + j) * m + k;
                    int left = ((i - 1) * n + j) * m + k;
                    int right = ((i + 1) * n + j) * m + k;
                    int bottom = (i * n + (j - 1)) * m + k;
                    int top = (i * n + (j + 1)) * m + k;
                    int back = (i * n + j) * m + (k - 1);
                    int front = (i * n + j) * m + (k + 1);

                    float s = this->s[center];
                    float sx0 = this->s[left];
                    float sx1 = this->s[right];
                    float sy0 = this->s[bottom];
                    float sy1 = this->s[top];
                    float sz0 = this->s[back];
                    float sz1 = this->s[front];
                    s = sx0 + sx1 + sy0 + sy1 + sz0 + sz1;
                    if (s == 0.0)
                        continue;

                    float div = u[right] - u[center] +
                        v[top] - v[center] +
                        w[front] - w[back];

                    if (particleRestDensity > 0.0 && compensateDrift) {
                        float k = 1.0;
                        float compression = particleDensity[center] - particleRestDensity;
                        if (compression > 0.0)
                            div -= k * compression;
                    }

                    float pressure = -div / s;
                    pressure *= overRelaxation;
                    p[center] += cp * pressure;

                    u[center] -= sx0 * pressure;
                    u[right] += sx1 * pressure;
                    v[center] -= sy0 * pressure;
                    v[top] += sy1 * pressure;
                    w[center] -= sz0 * pressure;
                    w[front] += sz1 * pressure;
                }
            }
        }
    }
}

void FLIPSimulation::step()
{
    currentStep++;
//    float dt = 0.1; // Time step
    float gravityValue = 9.81; // Acceleration due to gravity
//    Vector3 obstaclePos(0, 0, 0);
//    float obstacleRadius = 0.0; // Radius of obstacle
    int numIterations = 100; // Number of iterations for each step
    float overRelaxation = 1.0; // Over-relaxation parameter
    bool compensateDrift = true; // Whether to compensate for drift
    Vector3 gravity(0, 0, -gravityValue);
    float flipRatio = 0.9;

//    std::vector<Particle> copy = {};
//    for (auto& particle : particles)
//        if ((particle.position - obstaclePos).norm2() >= obstacleRadius * obstacleRadius)
//            copy.push_back(particle);
//    particles = copy;


    savedState = this->particles;
    integrateParticles(dt, gravity);
    this->useVelocityForCollisionDetection = true;
    handleCollisions();
    pushParticlesApart(numIterations);
    this->useVelocityForCollisionDetection = false;
    handleCollisions();
    updateParticleDensity();
    transferVelocities(true, flipRatio);
    solveIncompressibility(numIterations, dt, overRelaxation, compensateDrift);
    transferVelocities(false, flipRatio);
}

void FLIPSimulation::simulate() {

    while (true) {
        step();
    }
}

Matrix3<Vector3> FLIPSimulation::getVelocities(int newSizeX, int newSizeY, int newSizeZ)
{
    if (_cachedStep != currentStep) {
        _cachedStep = currentStep;
        Matrix3<Vector3> velocities(dimensions);
        Matrix3<float> amount(dimensions);

        for (auto& particle : this->particles) {
            velocities[particle.position] += particle.velocity;
            amount[particle.position] += 1;
        }
        for (size_t i = 0; i < velocities.size(); i++)
            velocities[i] /= amount[i];

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



















/*
Scene::Scene() {
    this->setupScene();
}
void Scene::setupScene() {
    gravity = -9.8f;                // Set the gravity value
    dt = 0.01f;                     // Time step size
    flipRatio = 0.95f;              // Ratio of FLIP velocities to grid velocities
    numPressureIters = 20;          // Number of pressure solver iterations
    numParticleIters = 5;           // Number of particle advection iterations
    frameNr = 0;                    // Frame number
    overRelaxation = 1.0f;          // Over-relaxation factor for pressure solver
    compensateDrift = true;         // Flag to compensate particle drift
    separateParticles = true;       // Flag to separate particles that are too close
    obstacleX = 0.5f;               // X-coordinate of the obstacle center
    obstacleY = 0.5f;               // Y-coordinate of the obstacle center
    obstacleZ = 0.5f;               // Z-coordinate of the obstacle center
    obstacleRadius = 0.2f;          // Radius of the obstacle
    paused = false;                 // Flag to pause the simulation
    showObstacle = true;            // Flag to show the obstacle
    obstacleVelX = 0.0f;            // Velocity of the obstacle in the X-direction
    obstacleVelY = 0.0f;            // Velocity of the obstacle in the Y-direction
    obstacleVelZ = 0.0f;            // Velocity of the obstacle in the Z-direction
    showParticles = true;           // Flag to show the particles
    showGrid = true;                // Flag to show the grid

    // Create a new instance of FlipFluid
    float density = 1000.0f;        // Fluid density
    float width = 1.0f;             // Width of the fluid domain
    float height = 1.0f;            // Height of the fluid domain
    float depth = 1.0f;             // Depth of the fluid domain
    float spacing = 1.0f / 64.0f;   // Spacing between grid cells
    float particleRadius = spacing; // Radius of the particles
    int maxParticles = 10000;       // Maximum number of particles
    fluid = new FLIPSimulation(density, width, height, depth, spacing, particleRadius, maxParticles);
}

void Scene::setObstacle(float x, float y, float z, bool reset) {
    obstacleX = x;
    obstacleY = y;
    obstacleZ = z;

    if (reset) {
        // Reset the simulation if requested
        fluid->reset();
        frameNr = 0;
    }
}

void Scene::simulate() {
    if (!paused) {
        // Simulate one time step
        fluid->simulate(dt, gravity, flipRatio, numPressureIters, numParticleIters, overRelaxation,
                        compensateDrift, separateParticles, obstacleX, obstacleY, obstacleZ, obstacleRadius);

        // Update the frame number
        frameNr++;
    }
}

void Scene::update() {
    // Main simulation loop
//    while (true) {
        // Update the scene and render

        // Simulate one time step
        simulate();

        // Render the scene

        // Update the scene variables based on user input or other events
//    }
}

FLIPSimulation::FLIPSimulation(float density, float width, float height, float depth, float spacing, float particleRadius, int maxParticles)
{
    this->grid = Matrix3<Vector3>(width, depth, height);
    this->spacing = spacing;
    this->nbParticles = maxParticles;
    this->reset();
}

void FLIPSimulation::reset() {
    // Reset the simulation state
    particles.clear();
    velocities.clear();
    grid.reset();
//    activeCells.clear();

    particles.resize(this->nbParticles);
    for (auto& p : particles)
        p = Vector3::random(grid.getDimensions());
    velocities.resize(this->nbParticles);
}

void FLIPSimulation::simulate(float dt, float gravity, float flipRatio, int numPressureIters,
                         int numParticleIters, float overRelaxation, bool compensateDrift,
                         bool separateParticles, float obstacleX, float obstacleY, float obstacleZ,
                         float obstacleRadius) {
    // Apply external forces
    applyGravity(gravity);

    // Transfer velocities from particles to the grid
    transferToGrid();

    // Apply velocity boundary conditions on the grid
    applyBoundaryConditions();

    // Solve pressure and update grid velocities
    solvePressure(numPressureIters, overRelaxation);

    // Update the particle velocities
    updateParticleVelocities(flipRatio);

    // Advect particles
    advectParticles(dt, numParticleIters, obstacleX, obstacleY, obstacleZ, obstacleRadius);

    // Separate particles that are too close
    if (separateParticles) {
        separateCloseParticles();
    }

    // Compensate for particle drift
    if (compensateDrift) {
        compensateParticleDrift();
    }
}

void FLIPSimulation::applyGravity(float gravity) {
    // Apply gravity to the velocities of the particles
    for (int i = 0; i < particles.size(); i++) {
        velocities[i].y += gravity * dt;
    }
}

void FLIPSimulation::transferToGrid() {
    // Clear the grid
    grid.clear();

    // Transfer velocities from particles to the grid cells
    for (int i = 0; i < particles.size(); i++) {
        Vector3 particle = particles[i];
        Vector3 velocity = velocities[i];

        // Determine the cell indices containing the particle
        int iCell = int(particle.x / spacing);
        int jCell = int(particle.y / spacing);
        int kCell = int(particle.z / spacing);

        // Compute the fractional distance within the cell
        float xFrac = (particle.x - iCell * spacing) / spacing;
        float yFrac = (particle.y - jCell * spacing) / spacing;
        float zFrac = (particle.z - kCell * spacing) / spacing;

        // Interpolate the velocity to the grid cells using trilinear interpolation
        grid(iCell, jCell, kCell) += (1.0f - xFrac) * (1.0f - yFrac) * (1.0f - zFrac) * velocity;
        grid(iCell + 1, jCell, kCell) += xFrac * (1.0f - yFrac) * (1.0f - zFrac) * velocity;
        grid(iCell, jCell + 1, kCell) += (1.0f - xFrac) * yFrac * (1.0f - zFrac) * velocity;
        grid(iCell, jCell, kCell + 1) += (1.0f - xFrac) * (1.0f - yFrac) * zFrac * velocity;
        grid(iCell + 1, jCell + 1, kCell) += xFrac * yFrac * (1.0f - zFrac) * velocity;
        grid(iCell + 1, jCell, kCell + 1) += xFrac * (1.0f - yFrac) * zFrac * velocity;
        grid(iCell, jCell + 1, kCell + 1) += (1.0f - xFrac) * yFrac * zFrac * velocity;
        grid(iCell + 1, jCell + 1, kCell + 1) += xFrac * yFrac * zFrac * velocity;
    }
}

void FLIPSimulation::applyBoundaryConditions() {
    // Apply boundary conditions to the velocities on the grid
    int numCellsX = grid.width();
    int numCellsY = grid.height();
    int numCellsZ = grid.depth();

    // Apply boundary conditions to the velocity components along the x-axis
    for (int j = 0; j < numCellsY; j++) {
        for (int k = 0; k < numCellsZ; k++) {
            // Left boundary
            grid(0, j, k).x = 0.0f;
            grid(0, j, k).y = grid(1, j, k).y;
            grid(0, j, k).z = grid(1, j, k).z;

            // Right boundary
            grid(numCellsX - 1, j, k).x = 0.0f;
            grid(numCellsX - 1, j, k).y = grid(numCellsX - 2, j, k).y;
            grid(numCellsX - 1, j, k).z = grid(numCellsX - 2, j, k).z;
        }
    }

    // Apply boundary conditions to the velocity components along the y-axis
    for (int i = 0; i < numCellsX; i++) {
        for (int k = 0; k < numCellsZ; k++) {
            // Bottom boundary
            grid(i, 0, k).x = grid(i, 1, k).x;
            grid(i, 0, k).y = 0.0f;
            grid(i, 0, k).z = grid(i, 1, k).z;

            // Top boundary
            grid(i, numCellsY - 1, k).x = grid(i, numCellsY - 2, k).x;
            grid(i, numCellsY - 1, k).y = 0.0f;
            grid(i, numCellsY - 1, k).z = grid(i, numCellsY - 2, k).z;
        }
    }

    // Apply boundary conditions to the velocity components along the z-axis
    for (int i = 0; i < numCellsX; i++) {
        for (int j = 0; j < numCellsY; j++) {
            // Back boundary
            grid(i, j, 0).x = grid(i, j, 1).x;
            grid(i, j, 0).y = grid(i, j, 1).y;
            grid(i, j, 0).z = 0.0f;

            // Front boundary
            grid(i, j, numCellsZ - 1).x = grid(i, j, numCellsZ - 2).x;
            grid(i, j, numCellsZ - 1).y = grid(i, j, numCellsZ - 2).y;
            grid(i, j, numCellsZ - 1).z = 0.0f;
        }
    }
}

void FLIPSimulation::solvePressure(int numIterations, float overRelaxation) {
    // Solve the pressure Poisson equation using the Gauss-Seidel method
    int numCellsX = grid.width();
    int numCellsY = grid.height();
    int numCellsZ = grid.depth();

    gridW = Matrix3<float>(numCellsX, numCellsY, numCellsZ);

    for (int iter = 0; iter < numIterations; iter++) {
        for (int i = 0; i < numCellsX; i++) {
            for (int j = 0; j < numCellsY; j++) {
                for (int k = 0; k < numCellsZ; k++) {
                    // Compute the index of the neighboring cells
                    int left = std::max(i - 1, 0);
                    int right = std::min(i + 1, numCellsX - 1);
                    int bottom = std::max(j - 1, 0);
                    int top = std::min(j + 1, numCellsY - 1);
                    int back = std::max(k - 1, 0);
                    int front = std::min(k + 1, numCellsZ - 1);

                    // Compute the divergence of the velocity field
                    float divergence =
                        (grid(right, j, k).x - grid(left, j, k).x) +
                        (grid(i, top, k).y - grid(i, bottom, k).y) +
                        (grid(i, j, front).z - grid(i, j, back).z);

                    // Compute the pressure correction using the Gauss-Seidel formula
                    float pressureCorrection =
                        (divergence - gridW(i, j, k)) / 6.0f * overRelaxation;

                    // Update the pressure and divergence values
                    gridW(i, j, k) += pressureCorrection;
                    grid(left, j, k).x -= pressureCorrection;
                    grid(right, j, k).x -= pressureCorrection;
                    grid(i, bottom, k).y -= pressureCorrection;
                    grid(i, top, k).y -= pressureCorrection;
                    grid(i, j, back).z -= pressureCorrection;
                    grid(i, j, front).z -= pressureCorrection;
                }
            }
        }
    }
}

void FLIPSimulation::updateParticleVelocities(float flipRatio) {
    // Update the velocities of the particles using the FLIP algorithm
    for (int i = 0; i < particles.size(); i++) {
        Vector3 particle = particles[i];

        // Determine the cell indices containing the particle
        int iCell = int(particle.x / spacing);
        int jCell = int(particle.y / spacing);
        int kCell = int(particle.z / spacing);

        // Compute the fractional distance within the cell
        float xFrac = (particle.x - iCell * spacing) / spacing;
        float yFrac = (particle.y - jCell * spacing) / spacing;
        float zFrac = (particle.z - kCell * spacing) / spacing;

        // Interpolate the velocity from the grid to the particle position
        Vector3 velocity = (1.0f - xFrac) * (1.0f - yFrac) * (1.0f - zFrac) * grid(iCell, jCell, kCell);
        velocity += xFrac * (1.0f - yFrac) * (1.0f - zFrac) * grid(iCell + 1, jCell, kCell);
        velocity += (1.0f - xFrac) * yFrac * (1.0f - zFrac) * grid(iCell, jCell + 1, kCell);
        velocity += (1.0f - xFrac) * (1.0f - yFrac) * zFrac * grid(iCell, jCell, kCell + 1);
        velocity += xFrac * yFrac * (1.0f - zFrac) * grid(iCell + 1, jCell + 1, kCell);
        velocity += xFrac * (1.0f - yFrac) * zFrac * grid(iCell + 1, jCell, kCell + 1);
        velocity += (1.0f - xFrac) * yFrac * zFrac * grid(iCell, jCell + 1, kCell + 1);
        velocity += xFrac * yFrac * zFrac * grid(iCell + 1, jCell + 1, kCell + 1);

        // Update the particle velocity using the FLIP ratio
        velocities[i] = flipRatio * velocity + (1.0f - flipRatio) * velocities[i];
    }
}

void FLIPSimulation::advectParticles(float dt, int numIterations, float obstacleX, float obstacleY,
                                float obstacleZ, float obstacleRadius) {
    // Advect the particles using the PIC/FLIP method
    float invSpacing = 1.0f / spacing;

    for (int iter = 0; iter < numIterations; iter++) {
        for (int i = 0; i < particles.size(); i++) {
            Vector3 particle = particles[i];
            Vector3 velocity = velocities[i];

            // Compute the cell indices containing the particle
            int iCell = int(particle.x * invSpacing);
            int jCell = int(particle.y * invSpacing);
            int kCell = int(particle.z * invSpacing);

            // Compute the fractional distance within the cell
            float xFrac = (particle.x - iCell * spacing) * invSpacing;
            float yFrac = (particle.y - jCell * spacing) * invSpacing;
            float zFrac = (particle.z - kCell * spacing) * invSpacing;

            // Interpolate the velocity from the grid to the particle position
            Vector3 gridVelocity =
                (1.0f - xFrac) * (1.0f - yFrac) * (1.0f - zFrac) * grid(iCell, jCell, kCell);
            gridVelocity += xFrac * (1.0f - yFrac) * (1.0f - zFrac) * grid(iCell + 1, jCell, kCell);
            gridVelocity += (1.0f - xFrac) * yFrac * (1.0f - zFrac) * grid(iCell, jCell + 1, kCell);
            gridVelocity += (1.0f - xFrac) * (1.0f - yFrac) * zFrac * grid(iCell, jCell, kCell + 1);
            gridVelocity += xFrac * yFrac * (1.0f - zFrac) * grid(iCell + 1, jCell + 1, kCell);
            gridVelocity += xFrac * (1.0f - yFrac) * zFrac * grid(iCell + 1, jCell, kCell + 1);
            gridVelocity += (1.0f - xFrac) * yFrac * zFrac * grid(iCell, jCell + 1, kCell + 1);
            gridVelocity += xFrac * yFrac * zFrac * grid(iCell + 1, jCell + 1, kCell + 1);

            // Update the particle position using advection
            particle += dt * gridVelocity;

            // Apply the obstacle constraints
            float distanceToObstacle = std::sqrt((particle.x - obstacleX) * (particle.x - obstacleX) +
                                                 (particle.y - obstacleY) * (particle.y - obstacleY) +
                                                 (particle.z - obstacleZ) * (particle.z - obstacleZ));
            if (distanceToObstacle < obstacleRadius) {
                // Move the particle to the surface of the obstacle
                particle = Vector3(obstacleX, obstacleY, obstacleZ) +
                           obstacleRadius * (particle - Vector3(obstacleX, obstacleY, obstacleZ)) /
                               distanceToObstacle;
            }

            // Update the particle position and velocity
            particles[i] = particle;
            velocities[i] = velocity;
        }
    }
}

void FLIPSimulation::separateCloseParticles() {
    // Separate particles that are too close to each other
    float separationDistance = 0.5f * spacing;

    for (int i = 0; i < particles.size(); i++) {
        Vector3 particle = particles[i];

        for (int j = i + 1; j < particles.size(); j++) {
            Vector3 otherParticle = particles[j];

            if ((particle - otherParticle).length() < separationDistance) {
                // Move the particles apart
                Vector3 separationDirection = (particle - otherParticle).normalized();
                particles[i] += 0.5f * separationDistance * separationDirection;
                particles[j] -= 0.5f * separationDistance * separationDirection;
            }
        }
    }
}

void FLIPSimulation::compensateParticleDrift() {
    // Compensate for the particle drift caused by advection
    for (int i = 0; i < particles.size(); i++) {
        Vector3 particle = particles[i];
        Vector3 velocity = velocities[i];

        // Compute the cell indices containing the particle
        int iCell = int(particle.x / spacing);
        int jCell = int(particle.y / spacing);
        int kCell = int(particle.z / spacing);

        // Compute the fractional distance within the cell
        float xFrac = (particle.x - iCell * spacing) / spacing;
        float yFrac = (particle.y - jCell * spacing) / spacing;
        float zFrac = (particle.z - kCell * spacing) / spacing;

        // Interpolate the velocity from the grid to the particle position
        Vector3 gridVelocity =
            (1.0f - xFrac) * (1.0f - yFrac) * (1.0f - zFrac) * grid(iCell, jCell, kCell);
        gridVelocity += xFrac * (1.0f - yFrac) * (1.0f - zFrac) * grid(iCell + 1, jCell, kCell);
        gridVelocity += (1.0f - xFrac) * yFrac * (1.0f - zFrac) * grid(iCell, jCell + 1, kCell);
        gridVelocity += (1.0f - xFrac) * (1.0f - yFrac) * zFrac * grid(iCell, jCell, kCell + 1);
        gridVelocity += xFrac * yFrac * (1.0f - zFrac) * grid(iCell + 1, jCell + 1, kCell);
        gridVelocity += xFrac * (1.0f - yFrac) * zFrac * grid(iCell + 1, jCell, kCell + 1);
        gridVelocity += (1.0f - xFrac) * yFrac * zFrac * grid(iCell, jCell + 1, kCell + 1);
        gridVelocity += xFrac * yFrac * zFrac * grid(iCell + 1, jCell + 1, kCell + 1);

        // Compute the difference between the particle velocity and the grid velocity
        Vector3 velocityDifference = velocity - gridVelocity;

        // Update the particle velocity to compensate for the drift
        velocities[i] -= velocityDifference;
    }
}
*/
