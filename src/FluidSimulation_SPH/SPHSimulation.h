#ifndef SPHSIMULATION_H
#define SPHSIMULATION_H

#include <vector>
#include "DataStructure/Matrix3.h"
#include "DataStructure/Vector3.h"

namespace SPHSimulation
{
// Define a structure for a particle in the system.
class Particle {
public:
    Particle(float x, float y, float z)
        : x(x), y(y), z(z)
    {}

    float x, y, z;  // Position.
    float vx, vy, vz;  // Velocity.
    float force_x, force_y, force_z;  // Acceleration.
    float density, pressure, mass;  // Physical properties.

    std::vector<Particle*> neighbors;
};

float distance(Particle& p1, Particle& p2);

// Define the SPH system.
class SPHSystem {
    private:
    std::vector<Particle> particles;  // The particles in the system.
    float h;  // The smoothing length.
    float k;  // Gas constant.
    float rho0;  // Reference density.
    float gamma;  // Polytopic constant for water

    public:
    SPHSystem(float h, float k, float rho0) : h(h), k(k), rho0(rho0) {
        k = 3e6;  // Bulk modulus for water
        rho0 = 1000;  // Rest density for water
        gamma = 7;  // Polytopic constant for water
    }

    // Add a particle to the system.
    void addParticle(Particle p) {
        particles.push_back(p);
    }

    float kernelFunction(float dist, float maxDist) {
        float ratio = dist / maxDist;
        return std::max(std::min(1.f - ratio * ratio, 1.f), 0.f);
    }
    float derivativeKernelFunction(float dist, float maxDist) {
        return -(2.0 * dist) / (maxDist * maxDist); // Function directly defined after the "kernelFunction" formula
    }

    // The main simulation loop.
    void simulate() {
        float dt = 0.01;  // Time step.

        // First, compute the density and pressure for all particles.
        for (int i = 0; i < particles.size(); i++) {
          computeDensity(particles[i]);
          computePressure(particles[i]);
        }

        // Then compute the forces for all particles.
        for (int i = 0; i < particles.size(); i++) {
          computeForces(particles[i]);
        }

        // Then, do the time integration using the Leapfrog-Verlet method.
        for (int i = 0; i < particles.size(); i++) {
          // Half step update for the velocity.
          particles[i].vx += particles[i].force_x * dt / 2;
          particles[i].vy += particles[i].force_y * dt / 2;
          particles[i].vz += particles[i].force_z * dt / 2;

          // Full step update for the position.
          particles[i].x += particles[i].vx * dt;
          particles[i].y += particles[i].vy * dt;
          particles[i].z += particles[i].vz * dt;

          // Re-compute the forces with the new position.
          computeForces(particles[i]);

          // Complete the velocity update with the new force.
          particles[i].vx += particles[i].force_x * dt / 2;
          particles[i].vy += particles[i].force_y * dt / 2;
          particles[i].vz += particles[i].force_z * dt / 2;
        }
    }

    // Compute the density for a particle.
    void computeDensity(Particle& p) {
      float density = 0;

      /*
       * If at one point I finish the spatial partitionning...
      float cellSize = 1.f;

      // Loop over all particles in the same cell and neighboring cells.
      for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
          for (int dz = -1; dz <= 1; dz++) {
            int x = p.x / cellSize + dx;
            int y = p.y / cellSize + dy;
            int z = p.z / cellSize + dz;
            if (grid.count(x) && grid[x].count(y) && grid[x][y].count(z)) {
              for (Particle* neighbor : grid[x][y][z].particles) {
                float r = distance(p, *neighbor);
                if (r < h) {
                  float W = kernelFunction(r, h);  // Kernel function depends on your SPH formulation.
                  density += neighbor->mass * W;
                }
              }
            }
          }
        }
      }
      */

      for (Particle& neighbor : particles) {
        float r = distance(p, neighbor);
        if (r < h) {
          float W = kernelFunction(r, h);  // Kernel function depends on your SPH formulation.
          density += neighbor.mass * W;
        }
      }

      p.density = density;
    }

    // Compute the pressure for a particle.
    void computePressure(Particle& p) {

      p.pressure = k * (pow(p.density / rho0, gamma) - 1);
    }

    // Compute the forces for a particle.
    void computeForces(Particle& p) {
      // ... (pressure force calculation as before)
        float fx, fy, fz;

      // Viscosity force
      float mu = 0.01;  // Viscosity coefficient
      for (Particle* neighbor : p.neighbors) {
        float r = distance(p, *neighbor);
        if (r < h) {
          float W = kernelFunction(r, h);
          float dw = derivativeKernelFunction(r, h);
          float q = r / h;
          float u = (p.vx - neighbor->vx) * (p.x - neighbor->x) + (p.vy - neighbor->vy) * (p.y - neighbor->y) + (p.vz - neighbor->vz) * (p.z - neighbor->z);
          if (u < 0) {
            float pi = -mu * (1 - q) * (u / (r + 0.01 * h)) * (W / (p.density * neighbor->density));
            fx += pi * (p.x - neighbor->x) / r;
            fy += pi * (p.y - neighbor->y) / r;
            fz += pi * (p.z - neighbor->z) / r;
          }
        }
      }

      p.force_x = fx;
      p.force_y = fy;
      p.force_z = fz;
    }

    void computeVelocityGradient(Particle& p) {
        float cellSize = 1.0;

      float dvx_dx = 0, dvx_dy = 0, dvx_dz = 0;
      float dvy_dx = 0, dvy_dy = 0, dvy_dz = 0;
      float dvz_dx = 0, dvz_dy = 0, dvz_dz = 0;

      /*
       * If at one point I finish the spatial partitionning...
      // Loop over all particles in the same cell and neighboring cells.
      for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
          for (int dz = -1; dz <= 1; dz++) {
            int x = p.x / cellSize + dx;
            int y = p.y / cellSize + dy;
            int z = p.z / cellSize + dz;
            if (grid.count(x) && grid[x].count(y) && grid[x][y].count(z)) {
              for (Particle* neighbor : grid[x][y][z].particles) {
                float r = distance(p, *neighbor);
                if (r < h) {
                  float W = kernelFunction(r, h);
                  float dw_dx = derivativeKernelFunction(r, h) * (p.x - neighbor->x) / r;
                  float dw_dy = derivativeKernelFunction(r, h) * (p.y - neighbor->y) / r;
                  float dw_dz = derivativeKernelFunction(r, h) * (p.z - neighbor->z) / r;

                  dvx_dx += neighbor->mass * (neighbor->vx - p.vx) * dw_dx / neighbor->density;
                  dvx_dy += neighbor->mass * (neighbor->vx - p.vx) * dw_dy / neighbor->density;
                  dvx_dz += neighbor->mass * (neighbor->vx - p.vx) * dw_dz / neighbor->density;

                  dvy_dx += neighbor->mass * (neighbor->vy - p.vy) * dw_dx / neighbor->density;
                  dvy_dy += neighbor->mass * (neighbor->vy - p.vy) * dw_dy / neighbor->density;
                  dvy_dz += neighbor->mass * (neighbor->vy - p.vy) * dw_dz / neighbor->density;

                  dvz_dx += neighbor->mass * (neighbor->vz - p.vz) * dw_dx / neighbor->density;
                  dvz_dy += neighbor->mass * (neighbor->vz - p.vz) * dw_dy / neighbor->density;
                  dvz_dz += neighbor->mass * (neighbor->vz - p.vz) * dw_dz / neighbor->density;
                }
              }
            }
          }
        }
      }
      */
      for (Particle& neighbor : particles) {
        float r = distance(p, neighbor);
        if (r < h) {
          float W = kernelFunction(r, h);
          float dw_dx = derivativeKernelFunction(r, h) * (p.x - neighbor.x) / r;
          float dw_dy = derivativeKernelFunction(r, h) * (p.y - neighbor.y) / r;
          float dw_dz = derivativeKernelFunction(r, h) * (p.z - neighbor.z) / r;

          dvx_dx += neighbor.mass * (neighbor.vx - p.vx) * dw_dx / neighbor.density;
          dvx_dy += neighbor.mass * (neighbor.vx - p.vx) * dw_dy / neighbor.density;
          dvx_dz += neighbor.mass * (neighbor.vx - p.vx) * dw_dz / neighbor.density;

          dvy_dx += neighbor.mass * (neighbor.vy - p.vy) * dw_dx / neighbor.density;
          dvy_dy += neighbor.mass * (neighbor.vy - p.vy) * dw_dy / neighbor.density;
          dvy_dz += neighbor.mass * (neighbor.vy - p.vy) * dw_dz / neighbor.density;

          dvz_dx += neighbor.mass * (neighbor.vz - p.vz) * dw_dx / neighbor.density;
          dvz_dy += neighbor.mass * (neighbor.vz - p.vz) * dw_dy / neighbor.density;
          dvz_dz += neighbor.mass * (neighbor.vz - p.vz) * dw_dz / neighbor.density;
        }
      }
    }

    // This function computes the shear stress vector given the velocity gradient and the surface normal.
    Vector3 computeShearStress(Matrix3<float>& velocityGradient, Vector3& surfaceNormal) {
      float dynamicViscosity = 0.001;  // Dynamic viscosity of the fluid, in Pa*s.

      // Compute the strain rate tensor, which is the symmetric part of the velocity gradient tensor.
      Matrix3<float> strainRate;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          strainRate(i, j) = 0.5 * (velocityGradient(i, j) + velocityGradient(j, i));
        }
      }

      // Compute the shear rate tensor, which is the strain rate tensor with the diagonal elements subtracted out.
      Matrix3 shearRate = strainRate;
      for (int i = 0; i < 3; i++) {
        shearRate(i, i) -= strainRate.trace() / 3.0;
      }

      // Multiply the shear rate tensor by the dynamic viscosity to get the shear stress tensor.
      Matrix3 shearStressTensor = shearRate * dynamicViscosity;
      Matrix3<float>& m = shearStressTensor;
      Vector3& n = surfaceNormal;

      // Compute the shear stress vector by multiplying the shear stress tensor by the surface normal vector.
      Vector3 shearStress = Vector3(
                  m(0, 0) * n.x + m(0, 1) * n.y + m(0, 2) * n.z,
                  m(1, 0) * n.x + m(1, 1) * n.y + m(1, 2) * n.z,
                  m(2, 0) * n.x + m(2, 1) * n.y + m(2, 2) * n.z
                  );


//              shearStressTensor * surfaceNormal;

      return shearStress;
    }
};

}

#endif // SPHSIMULATION_H
