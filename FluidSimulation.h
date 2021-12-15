#ifndef FLUIDSIMULATION_H
#define FLUIDSIMULATION_H

#include "Vector3.h"

class FluidSimulation
{
public:
    FluidSimulation();
    FluidSimulation(int sizeX, int sizeY, int sizeZ, float dt, float diffusionAmount, float viscosity, int iterations);

    void setObstacles(std::vector<std::vector<std::vector<bool>>> new_obstacles);
    void addDensity(int x, int y, int z, float amount);
    void addVelocity(int x, int y, int z, Vector3 amount);

    std::vector<std::vector<std::vector<Vector3>>> getVelocities(int rescaleX = -1, int rescaleY = -1, int rescaleZ = -1);

    void step();
    void velocityStep();
    void densityStep();
    void project();

    int sizeX, sizeY, sizeZ;
    int scaleX, scaleY, scaleZ;
    float dt;
    float diffusionAmount;
    float viscosity;

    int iterations;

    std::vector<std::vector<std::vector<float>>> density_old;
    std::vector<std::vector<std::vector<float>>> density;

    std::vector<std::vector<std::vector<Vector3>>> velocity;
    std::vector<std::vector<std::vector<Vector3>>> velocity_old;

    std::vector<std::vector<std::vector<float>>> divergence;
    std::vector<std::vector<std::vector<float>>> pressure;

    std::vector<std::vector<std::vector<bool>>> obstacles;


    template<typename T>
    void swapArrays(std::vector<std::vector<std::vector<T> > >& arr, std::vector<std::vector<std::vector<T> > >& arr2)
    {
        std::vector<std::vector<std::vector<T> > > tmpArray = arr;
        arr = arr2;
        arr2 = tmpArray;
    }

    template<typename T>
    void diffuse(std::vector<std::vector<std::vector<T> > >& arr, std::vector<std::vector<std::vector<T> > >& arr2, float diff, bool inverseOnBounds = false)
    {
        // TODO : check "a"
        float a = dt * diff * (arr.size() - 2) * (arr2.size() - 2);
        this->solve_linear(arr, arr2, a, 1 + 6 * a, inverseOnBounds);
    }
    template<typename T>
    void set_bounds(std::vector<std::vector<std::vector<T> > >& arr, bool inverseOnBounds = false, bool nullifyOnBounds = true)
    {
        std::vector<std::vector<std::vector<T> > > tmpArray = arr;
        for (int x = 1; x < this->sizeX - 1; x++) {
            for (int y = 1; y < this->sizeY - 1; y++) {
                for (int z = 1; z < this->sizeZ - 1; z++) {
                    if(this->obstacles[x][y][z]) {
                        bool freeCellFound = false;
                        for (int dx = -1; dx <= 1; dx++) {
                            for (int dy = -1; dy <= 1; dy++) {
                                for (int dz = -1; dz <= 1; dz++) {
                                    if (dx != 0 || dy != 0 || dz != 0) {
                                        if (!obstacles[x + dx][y + dy][z + dz]) {
                                            tmpArray[x][y][z] = arr[x + dx][y + dy][z + dz] * (inverseOnBounds ? -1.0 : 1.0);
                                            freeCellFound = true;
                                        }
                                    }
                                }
                            }
                        }
                        if (!freeCellFound)
                            tmpArray[x][y][z] = T();
                    }
                }
            }
        }
        arr = tmpArray;

        int xSize = arr.size(), ySize = arr[0].size(), zSize = arr[0][0].size();
        int X1 = xSize - 1, Y1 = ySize - 1, Z1 = zSize - 1;
        int X2 = xSize - 2, Y2 = ySize - 2, Z2 = zSize - 2;
        for (int x = 1; x < this->sizeX - 1; x++) {
            for (int y = 1; y < this->sizeY - 1; y++) {
                arr[x][y][0] = arr[x][y][1] * ((inverseOnBounds ? -1.0 : 1.0) * (nullifyOnBounds ? 0.0 : 1.0));
                arr[x][y][zSize - 1] = arr[x][y][zSize - 2] * ((inverseOnBounds ? -1.0 : 1.0) * (nullifyOnBounds ? 0.0 : 1.0));
            }
        }
        for (int x = 1; x < this->sizeX - 1; x++) {
            for (int z = 1; z < this->sizeZ - 1; z++) {
                arr[x][0][z] = arr[x][1][z] * ((inverseOnBounds ? -1.0 : 1.0) * (nullifyOnBounds ? 0.0 : 1.0));
                arr[x][ySize - 1][z] = arr[x][ySize - 2][z] * ((inverseOnBounds ? -1.0 : 1.0) * (nullifyOnBounds ? 0.0 : 1.0));
            }
        }
        for (int y = 1; y < this->sizeY - 1; y++) {
            for (int z = 1; z < this->sizeZ - 1; z++) {
                arr[0][y][z] = arr[1][y][z] * ((inverseOnBounds ? -1.0 : 1.0) * (nullifyOnBounds ? 0.0 : 1.0));
                arr[xSize - 1][y][z] = arr[xSize - 2][y][z] * ((inverseOnBounds ? -1.0 : 1.0) * (nullifyOnBounds ? 0.0 : 1.0));
            }
        }

        arr[0 ][0 ][0 ] = ( arr[0 ][0 ][1 ] +
                            arr[0 ][1 ][0 ] +
                            arr[1 ][0 ][0 ]) * (1/3.f);
        arr[0 ][0 ][Z1] = ( arr[0 ][0 ][Z2] +
                            arr[0 ][1 ][Z1] +
                            arr[1 ][0 ][Z1]) * (1/3.f);
        arr[0 ][Y1][0 ] = ( arr[0 ][Y1][1 ] +
                            arr[0 ][Y2][0 ] +
                            arr[1 ][Y1][0 ]) * (1/3.f);
        arr[0 ][Y1][Z1] = ( arr[0 ][Y1][Z2] +
                            arr[0 ][Y2][Z1] +
                            arr[1 ][Y1][Z1]) * (1/3.f);
        arr[X1][0 ][0 ] = ( arr[X1][0 ][1 ] +
                            arr[X1][1 ][0 ] +
                            arr[X2][0 ][0 ]) * (1/3.f);
        arr[X1][0 ][Z1] = ( arr[X1][0 ][Z2] +
                            arr[X1][1 ][Z1] +
                            arr[X2][0 ][Z1]) * (1/3.f);
        arr[X1][Y1][0 ] = ( arr[X1][Y1][1 ] +
                            arr[X1][Y2][0 ] +
                            arr[X2][Y1][0 ]) * (1/3.f);
        arr[X1][Y1][Z1] = ( arr[X1][Y1][Z2] +
                            arr[X1][Y2][Z1] +
                            arr[X2][Y1][Z1]) * (1/3.f);
    }

    template<typename T>
    void solve_linear(std::vector<std::vector<std::vector<T> > >& arr1, std::vector<std::vector<std::vector<T> > >& arr0, float a, float c, bool inverseOnBounds = false)
    {
        return;
        // TODO : check cRecip
        float cRecip = 1.0 / c;
        for (int k = 0; k < this->iterations; k++) {
            for (int x = 1; x < this->sizeX - 1; x++) {
                for (int y = 1; y < this->sizeY - 1; y++) {
                    for (int z = 1; z < this->sizeZ - 1; z++) {
                        arr1[x][y][z] = (arr0[x][y][z] +
                                ( arr1[x+1][y  ][z  ]
                                + arr1[x-1][y  ][z  ]
                                + arr1[x  ][y+1][z  ]
                                + arr1[x  ][y-1][z  ]
                                + arr1[x  ][y  ][z+1]
                                + arr1[x  ][y  ][z-1]) * a) * cRecip;
                    }
                }
            }
            this->set_bounds(arr1, inverseOnBounds);
        }
    }

    template<typename T>
    void advect(std::vector<std::vector<std::vector<T> > >& arr, std::vector<std::vector<std::vector<T> > >& old_array,
                std::vector<std::vector<std::vector<Vector3> > >& velocity_array, bool inverseOnBounds = false)
    {
        float x0, y0, z0, x1, y1, z1;

        float dtx = this->dt * (this->sizeX - 2);
        float dty = this->dt * (this->sizeY - 2);
        float dtz = this->dt * (this->sizeZ - 2);

        float s0, s1, t0, t1, u0, u1;
        Vector3 tmpVec;

        for (int x = 1; x < this->sizeX - 1; x++) {
            for (int y = 1; y < this->sizeY - 1; y++) {
                for (int z = 1; z < this->sizeZ - 1; z++) {
                    tmpVec = Vector3(x, y, z) - (velocity_array[x][y][z] * Vector3(dtx, dty, dtz));

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

                    arr[x][y][z] = ((old_array[x0i][y0i][z0i] * u0 + old_array[x0i][y0i][z1i] * u1) * t0 +
                                    (old_array[x0i][y1i][z0i] * u0 + old_array[x0i][y1i][z1i] * u1) * t1) * s0 +
                                   ((old_array[x1i][y0i][z0i] * u0 + old_array[x1i][y0i][z1i] * u1) * t0 +
                                    (old_array[x1i][y1i][z0i] * u0 + old_array[x1i][y1i][z1i] * u1) * t1) * s1;
//                    arr[x][y][z] = s0 * ( t0 * (u0 * arr2[x0i][y0i][z0i] + u1 * arr2[x0i][y0i][z1i])
//                                        + t1 * (u0 * arr2[x0i][y1i][z0i] + u1 * arr2[x0i][y1i][z1i]))
//                                 + s1 * ( t0 * (u0 * arr2[x1i][y0i][z0i] + u1 * arr2[x1i][y0i][z1i])
//                                        + t1 * (u0 * arr2[x1i][y1i][z0i] + u1 * arr2[x1i][y1i][z1i]));
                }
            }
        }
        this->set_bounds(arr, inverseOnBounds);
    }
};

#endif // FLUIDSIMULATION_H
