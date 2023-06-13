#ifndef FLIPSIMULATION_H
#define FLIPSIMULATION_H

#include <iostream>
#include <cmath>
#include <vector>
#include "DataStructure/Matrix3.h"
#include "DataStructure/Vector3.h"

namespace FLIP {

class FlipFluid {
private:
//    int fNumX;
//    int fNumY;
//    int fNumZ;
//    float h;
//    float fInvSpacing;
//    int fNumCells;
//    std::vector<float> u;
//    std::vector<float> v;
//    std::vector<float> w;
//    std::vector<float> du;
//    std::vector<float> dv;
//    std::vector<float> dw;
//    std::vector<float> prevU;
//    std::vector<float> prevV;
//    std::vector<float> prevW;
//    std::vector<float> p;
//    std::vector<float> s;
//    std::vector<int> cellType;
//    std::vector<float> cellColor;



    // Other variables and methods...
    std::vector<Vector3> particles;
    std::vector<Vector3> velocities;
    float particleRadius;
    float density;
    int maxParticles;

    Matrix3<Vector3> grid;

    float spacing = 1.f;

public:
    FlipFluid(float density, float width, float height, float depth, float spacing, float particleRadius, int maxParticles) {
        this->particles = std::vector<Vector3>(maxParticles);
        this->velocities = std::vector<Vector3>(maxParticles);
        this->particleRadius = particleRadius;
        this->density = density;
        reset();
    }

    void reset() {
        this->particles = std::vector<Vector3>(maxParticles);
        this->velocities = std::vector<Vector3>(maxParticles);
        this->grid.reset();
    }

    void simulate(float dt, float gravity, float flipRatio, int numPressureIters, int numParticleIters,
                  float overRelaxation, bool compensateDrift, bool separateParticles,
                  float obstacleX, float obstacleY, float obstacleZ, float obstacleRadius) {
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

    void applyGravity(float gravity) {
        // Apply gravity to the velocities of the particles
        for (int i = 0; i < particles.size(); i++) {
            velocities[i].y += gravity;
        }
    }

    void transferToGrid() {
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

    void applyBoundaryConditions() {
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

    void solvePressure(int numIterations, float overRelaxation) {
        // Solve the pressure Poisson equation using the Gauss-Seidel method
        int numCellsX = grid.width();
        int numCellsY = grid.height();
        int numCellsZ = grid.depth();

        Matrix3<float> pressureValue(grid.getDimensions(), 0.f);

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
                            (divergence - pressureValue(i, j, k)) / 6.0f * overRelaxation;

                        // Update the pressure and divergence values
                        pressureValue(i, j, k) += pressureCorrection;
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

    void updateParticleVelocities(float flipRatio) {
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

    void advectParticles(float dt, int numIterations, float obstacleX, float obstacleY,
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
                Vector3 gridVelocity = (1.0f - xFrac) * (1.0f - yFrac) * (1.0f - zFrac) * grid(iCell, jCell, kCell);
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

    void separateCloseParticles() {
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

    void compensateParticleDrift() {
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
    FlipFluid* fluid;

    Scene() {
        // Constructor implementation...
    }

    void setupScene() {
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
        fluid = new FlipFluid(density, width, height, depth, spacing, particleRadius, maxParticles);
    }

    void setObstacle(float x, float y, float z, bool reset) {
        obstacleX = x;
        obstacleY = y;
        obstacleZ = z;

        if (reset) {
            // Reset the simulation if requested
            fluid->reset();
            frameNr = 0;
        }
    }

    void simulate() {
        if (!paused) {
            // Simulate one time step
            fluid->simulate(dt, gravity, flipRatio, numPressureIters, numParticleIters, overRelaxation,
                            compensateDrift, separateParticles, obstacleX, obstacleY, obstacleZ, obstacleRadius);

            // Update the frame number
            frameNr++;
        }
    }

    void update() {
        // Main simulation loop
        while (true) {
            // Update the scene and render

            // Simulate one time step
            simulate();

            // Render the scene

            // Update the scene variables based on user input or other events
        }
    }
};

}
#endif // FLIPSIMULATION_H
