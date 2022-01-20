#ifndef FLUIDSIMULATION_H
#define FLUIDSIMULATION_H

#include "DataStructure/Vector3.h"
#include "DataStructure/Matrix3.h"

// The process is only partially understood in my head, see Josh Stam's paper for more understanding (Real-Time Fluid Dynamics for Games, 2003)
class FluidSimulation
{
public:
    FluidSimulation();
    FluidSimulation(int sizeX, int sizeY, int sizeZ, float dt, float diffusionAmount, float viscosity, int iterations);

    void setObstacles(Matrix3<float> new_obstacles);
    void addDensity(int x, int y, int z, float amount);
    void addVelocity(int x, int y, int z, Vector3 amount);
    void setMaxSpeed(float speed);

    Matrix3<Vector3> getVelocities(int rescaleX = -1, int rescaleY = -1, int rescaleZ = -1);

    void step();
    void diffuseVelocity();
    void advectVelocity();
    void velocityStep();
    void densityStep();
    void projectVelocity();

    void setVelocityBounds();

    int sizeX, sizeY, sizeZ;
    float dt;
    float diffusionAmount;
    float viscosity;
    float maxSpeed = 1.0;
    float maxSpeedSquared;

    int currentStep = 0;

    int iterations;

    Matrix3<float> density_old;
    Matrix3<float> density;

    Matrix3<Vector3> velocity;
    Matrix3<Vector3> velocity_old;

    Matrix3<float> divergence;
    Matrix3<float> pressure;

    Matrix3<float> obstacles;

    Vector3 meanVel;


    template<typename T>
    void swapArrays(Matrix3<T>& arr, Matrix3<T>& arr2)
    {
        Matrix3<T> tmpArray = arr;
        arr = arr2;
        arr2 = tmpArray;
    }

    template<typename T>
    void diffuse(Matrix3<T>& arr, Matrix3<T>& arr2, float diff, bool inverseOnBounds = false)
    {
        // TODO : check "a"
        float a = dt * diff * (arr.sizeX - 2) * (arr.sizeY - 2) * (arr.sizeZ - 2);
        this->solve_linear(arr, arr2, a, inverseOnBounds);
    }
    template<typename T>
    void set_bounds(Matrix3<T>& arr, bool inverseOnBounds = false, bool nullifyOnBounds = true)
    {
        nullifyOnBounds = true;
        inverseOnBounds = false;
        Matrix3<T> tmpArray = arr;

        for (int x = 1; x < this->sizeX - 1; x++) {
            for (int y = 1; y < this->sizeY - 1; y++) {
                for (int z = 1; z < this->sizeZ - 1; z++) {
                    if(this->obstacles(x, y, z) > .5) {
                        bool freeCellFound = false;
                        for (int dx = -1; dx <= 1; dx++) {
                            for (int dy = -1; dy <= 1; dy++) {
                                for (int dz = -1; dz <= 1; dz++) {
                                    if (dx != 0 || dy != 0 || dz != 0) {
                                        if (obstacles(x + dx, y + dy, z + dz) < .5) {
                                            tmpArray(x, y, z) = arr(x + dx, y + dy, z + dz) * (inverseOnBounds ? -1.0 : 1.0);
                                            freeCellFound = true;
                                        }
                                    }
                                }
                            }
                        }
                        if (!freeCellFound)
                            tmpArray(x, y, z) = T();
                    }
                }
            }
        }
        arr = tmpArray;

        int xSize = arr.sizeX, ySize = arr.sizeY, zSize = arr.sizeZ;
        int X1 = xSize - 1, Y1 = ySize - 1, Z1 = zSize - 1;
        int X2 = xSize - 2, Y2 = ySize - 2, Z2 = zSize - 2;
        for (int x = 1; x < this->sizeX - 1; x++) {
            for (int y = 1; y < this->sizeY - 1; y++) {
                arr(x, y, 0) = arr(x, y, 1) * ((inverseOnBounds ? -1.0 : 1.0) * (nullifyOnBounds ? 0.0 : 1.0));
                arr(x, y, zSize - 1) = arr(x, y, zSize - 2) * ((inverseOnBounds ? -1.0 : 1.0) * (nullifyOnBounds ? 0.0 : 1.0));
            }
        }
        for (int x = 1; x < this->sizeX - 1; x++) {
            for (int z = 1; z < this->sizeZ - 1; z++) {
                arr(x, 0, z) = arr(x, 1, z) * ((inverseOnBounds ? -1.0 : 1.0) * (nullifyOnBounds ? 0.0 : 1.0));
                arr(x, ySize - 1, z) = arr(x, ySize - 2, z) * ((inverseOnBounds ? -1.0 : 1.0) * (nullifyOnBounds ? 0.0 : 1.0));
            }
        }
        for (int y = 1; y < this->sizeY - 1; y++) {
            for (int z = 1; z < this->sizeZ - 1; z++) {
                arr(0, y, z) = arr(1, y, z) * ((inverseOnBounds ? -1.0 : 1.0) * (nullifyOnBounds ? 0.0 : 1.0));
                arr(xSize - 1, y, z) = arr(xSize - 2, y, z) * ((inverseOnBounds ? -1.0 : 1.0) * (nullifyOnBounds ? 0.0 : 1.0));
            }
        }

        arr(0 , 0 , 0 ) = ( arr(0 , 0 , 1 ) +
                            arr(0 , 1 , 0 ) +
                            arr(1 , 0 , 0 )) * (1/3.f);
        arr(0 , 0 , Z1) = ( arr(0 , 0 , Z2) +
                            arr(0 , 1 , Z1) +
                            arr(1 , 0 , Z1)) * (1/3.f);
        arr(0 , Y1, 0 ) = ( arr(0 , Y1, 1 ) +
                            arr(0 , Y2, 0 ) +
                            arr(1 , Y1, 0 )) * (1/3.f);
        arr(0 , Y1, Z1) = ( arr(0 , Y1, Z2) +
                            arr(0 , Y2, Z1) +
                            arr(1 , Y1, Z1)) * (1/3.f);
        arr(X1, 0 , 0 ) = ( arr(X1, 0 , 1 ) +
                            arr(X1, 1 , 0 ) +
                            arr(X2, 0 , 0 )) * (1/3.f);
        arr(X1, 0 , Z1) = ( arr(X1, 0 , Z2) +
                            arr(X1, 1 , Z1) +
                            arr(X2, 0 , Z1)) * (1/3.f);
        arr(X1, Y1, 0 ) = ( arr(X1, Y1, 1 ) +
                            arr(X1, Y2, 0 ) +
                            arr(X2, Y1, 0 )) * (1/3.f);
        arr(X1, Y1, Z1) = ( arr(X1, Y1, Z2) +
                            arr(X1, Y2, Z1) +
                            arr(X2, Y1, Z1)) * (1/3.f);
    }

    template<typename T>
    void solve_linear(Matrix3<T>& arr1, Matrix3<T>& arr0, float a, bool inverseOnBounds = false)
    {
        for (int k = 0; k < this->iterations; k++) {
            for (int x = 1; x < this->sizeX - 1; x++) {
                for (int y = 1; y < this->sizeY - 1; y++) {
                    for (int z = 1; z < this->sizeZ - 1; z++) {
                        arr1(x, y, z) = (arr0(x, y, z) +
                                ( arr1(x+1, y  , z  )
                                + arr1(x-1, y  , z  )
                                + arr1(x  , y+1, z  )
                                + arr1(x  , y-1, z  )
                                + arr1(x  , y  , z+1)
                                + arr1(x  , y  , z-1)) * a) / (1 + 6*a);
                    }
                }
            }
            this->set_bounds(arr1, inverseOnBounds);
        }
    }

    template<typename T>
    void advect(Matrix3<T>& arr, Matrix3<T>& old_array,
                Matrix3<Vector3>& velocity_array, bool inverseOnBounds = false)
    {
        float x0, y0, z0, x1, y1, z1;

//        float dtx = this->dt * (this->sizeX - 2);
//        float dty = this->dt * (this->sizeY - 2);
//        float dtz = this->dt * (this->sizeZ - 2);
        Vector3 dtVector = Vector3(this->sizeX, this->sizeY, this->sizeZ) * this->dt;

        float s0, s1, t0, t1, u0, u1;
        Vector3 tmpVec;

        for (int x = 1; x < this->sizeX - 1; x++) {
            for (int y = 1; y < this->sizeY - 1; y++) {
                for (int z = 1; z < this->sizeZ - 1; z++) {
                    tmpVec = Vector3(x, y, z) - (velocity_array(x, y, z) * dtVector); // Vector3(dtx, dty, dtz));

                    tmpVec.x = std::min(std::max(tmpVec.x, .5f), float(this->sizeX) + .5f);
                    x0 = std::floor(tmpVec.x);
                    x1 = x0 + 1.f;
                    s1 = tmpVec.x - x0;
                    s0 = 1.0 - s1;

                    tmpVec.y = std::min(std::max(tmpVec.y, .5f), float(this->sizeY) + .5f);
                    y0 = std::floor(tmpVec.y);
                    y1 = y0 + 1.f;
                    t1 = tmpVec.y - y0;
                    t0 = 1.0 - t1;

                    tmpVec.z = std::min(std::max(tmpVec.z, .5f), float(this->sizeZ) + .5f);
                    z0 = std::floor(tmpVec.z);
                    z1 = z0 + 1.f;
                    u1 = tmpVec.z - z0;
                    u0 = 1.0 - u1;

                    // TODO : Correct the simplifications caused by high velocities in "velocity_array"
                    int x0i = std::min(x0, float(this->sizeX - 1));
                    int x1i = std::min(x1, float(this->sizeX - 1));
                    int y0i = std::min(y0, float(this->sizeY - 1));
                    int y1i = std::min(y1, float(this->sizeY - 1));
                    int z0i = std::min(z0, float(this->sizeZ - 1));
                    int z1i = std::min(z1, float(this->sizeZ - 1));

                    arr(x, y, z) = ((old_array(x0i, y0i, z0i) * u0 + old_array(x0i, y0i, z1i) * u1) * t0 +
                                    (old_array(x0i, y1i, z0i) * u0 + old_array(x0i, y1i, z1i) * u1) * t1) * s0 +
                                   ((old_array(x1i, y0i, z0i) * u0 + old_array(x1i, y0i, z1i) * u1) * t0 +
                                    (old_array(x1i, y1i, z0i) * u0 + old_array(x1i, y1i, z1i) * u1) * t1) * s1;
                }
            }
        }
        this->set_bounds(arr, inverseOnBounds);
    }
};

#endif // FLUIDSIMULATION_H
