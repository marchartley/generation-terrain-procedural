#ifndef FLIPSIMULATION_H
#define FLIPSIMULATION_H

#include <iostream>
#include <cmath>
#include <vector>
#include "DataStructure/Vector3.h"
#include "DataStructure/Matrix3.h"
#include "FluidSimulation/FluidSimulation.h"

namespace FLIP {

class Particle {
public:
    Vector3 position;
    Vector3 velocity;
};


class FLIPSimulation : public FluidSimulation
{
public:
    int fNumX;
    int fNumY;
    int fNumZ;
    float h;
    float fInvSpacing;
    int fNumCells;
    std::vector<float> u;
    std::vector<float> v;
    std::vector<float> w;
    std::vector<float> du;
    std::vector<float> dv;
    std::vector<float> dw;
    std::vector<float> prevU;
    std::vector<float> prevV;
    std::vector<float> prevW;
    std::vector<float> p;
    std::vector<float> s;
    std::vector<int> cellType;
    std::vector<float> cellColor;

    std::vector<Particle> particles;
    std::vector<int> numCellParticles;
    float density;
    int maxParticles;

    std::vector<float> particleDensity;
    float particleRestDensity;
    float particleRadius;

    float pInvSpacing;
    int pNumX;
    int pNumY;
    int pNumZ;
    int pNumCells;

    std::vector<int> firstCellParticle;
    std::vector<int> cellParticleIds;

    Vector3 dimensions;

    int SOLID_CELL = 0;
    int FLUID_CELL = 1;
    int AIR_CELL = 2;

    float dt;

    FLIPSimulation();
    FLIPSimulation(float density, float width, float depth, float height, float spacing, float particleRadius, float maxParticles, float dt);
    virtual ~FLIPSimulation() {}
    void init(float density, float width, float depth, float height, float spacing, float particleRadius, float maxParticles, float dt);

    void integrateParticles(double dt, Vector3 gravity);
    void pushParticlesApart(int numIters);
//    void handleParticleCollisions(Vector3 obstaclePos, float obstacleRadius);
    void updateParticleDensity();
    void transferVelocities(bool toGrid, float flipRatio);
    void solveIncompressibility(int numIters, float dt, float overRelaxation, bool compensateDrift = true);
    void step();
    void simulate();

    void handleCollisions();


    Matrix3<Vector3> getVelocities(int newSizeX, int newSizeY, int newSizeZ);
    void addVelocity(int x, int y, int z, Vector3 amount);

    bool useVelocityForCollisionDetection;

    std::vector<Particle> savedState;
};











/*
class FLIPSimulation {
public:

    // Other variables and methods...

    std::vector<Vector3> particles;
    std::vector<Vector3> velocities;
    Matrix3<Vector3> grid;
    Matrix3<float> gridW;
//    Matrix3<Vector3> activeCells;
    float spacing = 1.f;
    int nbParticles;
    float dt = 0.01f;

public:
    FLIPSimulation(float density, float width, float height, float depth, float spacing, float particleRadius, int maxParticles);
    void reset();
    void simulate(float dt, float gravity, float flipRatio, int numPressureIters,
                  int numParticleIters, float overRelaxation, bool compensateDrift,
                  bool separateParticles, float obstacleX, float obstacleY, float obstacleZ,
                  float obstacleRadius);
    void applyGravity(float gravity);
    void transferToGrid();
    void applyBoundaryConditions();
    void solvePressure(int numIterations, float overRelaxation);
    void updateParticleVelocities(float flipRatio);
    void advectParticles(float dt, int numIterations, float obstacleX, float obstacleY,
                         float obstacleZ, float obstacleRadius);
    void separateCloseParticles();
    void compensateParticleDrift();
};

class Scene {
public:
    float gravity;
    float dt;
    float flipRatio;
    int numPressureIters;
    int numParticleIters;
    int frameNr;
    float overRelaxation;
    bool compensateDrift;
    bool separateParticles;
    float obstacleX;
    float obstacleY;
    float obstacleZ;
    float obstacleRadius;
    bool paused;
    bool showObstacle;
    float obstacleVelX;
    float obstacleVelY;
    float obstacleVelZ;
    bool showParticles;
    bool showGrid;
    FLIPSimulation* fluid;

    Scene();
    void setupScene();
    void setObstacle(float x, float y, float z, bool reset);
    void simulate();
    void update();
};

*/




/*
#include <iostream>
#include <cmath>
#include <vector>
#include "DataStructure/Matrix3.h"
#include "DataStructure/Vector3.h"
#include "Graphics/CubeMesh.h"
#include <vector>

namespace FLIP {

// Define a Cell structure with pressure
struct Cell {
    float u = 0.f, v = 0.f, w = 0.f;  // Velocity
    float p = 0.f;  // Pressure
    // ...
};
// Define a Particle
struct Particle {
    float u = 0.f, v = 0.f, w = 0.f;  // Velocity
    float p = 0.f;  // Pressure
    float u_prev = 0.f, v_prev = 0.f, w_prev = 0.f;  // Velocity
    float x = 0.f, y = 0.f, z = 0.f;  // Position
    // ...
};

class FLIPSimulation {
public:
    FLIPSimulation() {
        // Constructor code here...
    }

    // Grid of cells
    std::vector<std::vector<std::vector<Cell>>> grid;

    // List of particles
    std::vector<Particle> particles;
    // Constants for pressure solver
    const int pressureIterations = 50;
    const float pressureEpsilon = 0.01f;
    const Vector3 gravity = Vector3(0, 0, -.98f);
    float dt = 0.01f;

    int gridSize = 10;

    // List of triangles
    std::vector<std::vector<Vector3>> triangles;

    // Function to apply forces
    void applyForces() {
        for (Particle &p : particles) {
            // Apply gravitational force
            Vector3 gravityForce = gravity;
            p.u += gravityForce.x;
            p.v += gravityForce.y;
            p.w += gravityForce.z;
        }
    }
    // Function to advect particles
    void advect() {
        for (Particle &p : particles) {
            // Update particle position based on velocity
            p.x += p.u * dt;
            p.y += p.v * dt;
            p.z += p.w * dt;

            // Leapfrog integration
            p.u += 0.5f * (p.u_prev + gravity.x) * dt;
            p.v += 0.5f * (p.v_prev + gravity.y) * dt;
            p.w += 0.5f * (p.w_prev + gravity.z) * dt;

            // Store current velocity for next time step
            p.u_prev = p.u;
            p.v_prev = p.v;
            p.w_prev = p.w;
        }
    }

    // Function to compute the divergence at a cell
    float divergence(int i, int j, int k) {
        // Central difference approximation of divergence
        float div = 0;
        if (i > 0) div += grid[i][j][k].u - grid[i-1][j][k].u;
        if (j > 0) div += grid[i][j][k].v - grid[i][j-1][k].v;
        if (k > 0) div += grid[i][j][k].w - grid[i][j][k-1].w;
        return div;
    }

    // Function to solve the pressure Poisson equation
    void solvePressure() {
        for (int iter = 0; iter < pressureIterations; iter++) {
            // Jacobi iteration
            for (int i = 0; i < gridSize; i++) {
                for (int j = 0; j < gridSize; j++) {
                    for (int k = 0; k < gridSize; k++) {
                        float div = divergence(i, j, k);
                        float pressure = grid[i][j][k].p;
                        if (i > 0) pressure += grid[i-1][j][k].p;
                        if (j > 0) pressure += grid[i][j-1][k].p;
                        if (k > 0) pressure += grid[i][j][k-1].p;
                        pressure /= 3;
                        pressure -= div * pressureEpsilon;
                        grid[i][j][k].p = pressure;
                    }
                }
            }
        }
    }

    // Function to subtract pressure gradient from velocity
    void subtractPressureGradient() {
        for (int i = 0; i < gridSize; i++) {
            for (int j = 0; j < gridSize; j++) {
                for (int k = 0; k < gridSize; k++) {
                    if (i > 0) grid[i][j][k].u -= (grid[i][j][k].p - grid[i-1][j][k].p) / pressureEpsilon;
                    if (j > 0) grid[i][j][k].v -= (grid[i][j][k].p - grid[i][j-1][k].p) / pressureEpsilon;
                    if (k > 0) grid[i][j][k].w -= (grid[i][j][k].p - grid[i][j][k-1].p) / pressureEpsilon;
                }
            }
        }
    }

    // Function to solve fluid equations
    void solveFluidEquations() {
        // Compute forces (like gravity)
         applyForces();

        // Advect velocities
         advect();

        // Solve the pressure equations
        solvePressure();

        // Subtract the pressure gradient to enforce incompressibility
        subtractPressureGradient();
    }

    // Function to handle collisions
    void handleCollisions() {
        for (Particle &p : particles) {
            Vector3 startPos(p.x, p.y, p.z);
            Vector3 endPos = startPos + Vector3(p.u, p.v, p.w) * dt; // Assuming the velocity vector represents displacement for one time step.

            for (std::vector<Vector3> &tri : triangles) {
                Vector3 intersection = Collision::segmentToTriangleCollision(startPos, endPos, tri[0], tri[1], tri[2]);
                if (intersection.isValid()) {
                    // Collision detected, calculate the triangle normal
                    Vector3 edge1 = tri[1] - tri[0];
                    Vector3 edge2 = tri[2] - tri[0];
                    Vector3 normal = edge1.cross(edge2).normalized();

                    // Compute reflection vector
                    Vector3 velocity(p.u, p.v, p.w);
                    Vector3 reflection = velocity - 2 * velocity.dot(normal) * normal;

                    // Update particle velocity
                    p.u = reflection.x;
                    p.v = reflection.y;
                    p.w = reflection.z;

                    p.x = intersection.x + (intersection - endPos).x;
                    p.y = intersection.y + (intersection - endPos).y;
                    p.z = intersection.z + (intersection - endPos).z;
                }
            }
        }
    }

    // Main simulation loop
    void simulate() {
        while (true) {
            // Update fluid state
            solveFluidEquations();

            // Handle collisions
            handleCollisions();

            // Render or output the fluid state
            // ...

            // Advance to the next time step
            // ...
        }
    }

    void step() {
        solveFluidEquations();
        handleCollisions();
    }

    void initialize() {
        int nbParticles = 1000;
        this->particles.resize(nbParticles);
        for (int i = 0; i < nbParticles; i++) {
            particles[i] = Particle();
            particles[i].x = random_gen::generate(gridSize / 2, gridSize);
            particles[i].y = random_gen::generate(gridSize);
            particles[i].z = random_gen::generate(gridSize);
        }

        grid = std::vector<std::vector<std::vector<Cell>>>(gridSize, std::vector<std::vector<Cell>>(gridSize, std::vector<Cell>(gridSize)));

        Mesh m;
        std::vector<Vector3> verts = CubeMesh::cubesVertices;
        for (auto& v : verts)
            v = v * (gridSize - 2.f) + Vector3(1, 1, 1);
        m.fromArray(verts);
        triangles = m.getTriangles();
    }
};
*/
} // namespace FLIP
#endif // FLIPSIMULATION_H
