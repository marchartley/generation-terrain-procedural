#include "math.h"

void linearSolve(FlowGrid &x, const FlowGrid &x_0, Scalar a, Scalar c, const Indices &dim,
                 BoundarySetter setBoundaries, unsigned int iterations) {
    x = x_0;
    if (a == 0) {
        setBoundaries(x);
        return;
    }

    FlowGrid temp = x_0;
    for (unsigned int iteration = 0; iteration < iterations; ++iteration) {
#pragma omp parallel for
        for (FlowGrid::Index i = 1; i <= dim(0); ++i) {
            for (FlowGrid::Index j = 1; j <= dim(1); ++j) {
                for (FlowGrid::Index k = 1; k <= dim(2); ++k) {
                    temp(i, j, k) = (x_0(i, j, k) +
                                     a * (x(i - 1, j, k) + x(i + 1, j, k) +
                                          x(i, j - 1, k) + x(i, j + 1, k) +
                                          x(i, j, k - 1) + x(i, j, k + 1))) / c;
                }
            }
        }
        x = temp;

        setBoundaries(x);
   }
}

Scalar interpolate(const FlowGrid &grid, Location x) {
    Indices i = x.cast<FlowGrid::Index>();
    Indices j = i + 1;
    Location t = x - i.cast<Scalar>();
    Location s = 1 - t;
    return (s[2] * (s[1] * (s[0] * grid(i[0], i[1], i[2]) +
                            t[0] * grid(j[0], i[1], i[2])) +
                    t[1] * (s[0] * grid(i[0], j[1], i[2]) +
                            t[0] * grid(j[0], j[1], i[2]))) +
            t[2] * (s[1] * (s[0] * grid(i[0], i[1], j[2]) +
                            t[0] * grid(j[0], i[1], j[2])) +
                    t[1] * (s[0] * grid(i[0], j[1], j[2]) +
                            t[0] * grid(j[0], j[1], j[2]))));
}

void setBoundaries(FlowGrid &grid, int b, const Indices &dim) {
    for(FlowGrid::Index j = 1; j <= dim(1); ++j) {
        for (FlowGrid::Index k = 1; k <= dim(2); ++k) {
            grid(0, j, k) = (b == 0 ? -1 : 1) * grid(1, j, k);
            grid(dim(0) + 1, j, k) = (b == 0 ? -1 : 1) * grid(dim(0), j, k);
        }
    }
    for (FlowGrid::Index i = 1; i <= dim(0); ++i) {
        for (FlowGrid::Index k = 1; k <= dim(2); ++k) {
            grid(i, 0, k) = (b == 1 ? -1 : 1) * grid(i, 1, k);
            grid(i, dim(1) + 1, k) = (b == 1 ? -1 : 1) * grid(i, dim(1), k);
        }
    }
    for (FlowGrid::Index i = 1; i <= dim(0); ++i) {
        for (FlowGrid::Index j = 1; j <= dim(1); ++j) {
            grid(i, j, 0) = (b == 2 ? -1 : 1) * grid(i, j, 1);
            grid(i, j, dim(2) + 1) = (b == 2 ? -1 : 1) * grid(i, j, dim(2));
        }
    }

    grid(0, 0, 0) = (grid(1, 0, 0) + grid(0, 1, 0) + grid(0, 0, 1)) / 3;
    grid(dim(0) + 1, 0, 0) = (grid(dim(0), 0, 0) + grid(dim(0) + 1, 1, 0) +
                              grid(dim(0) + 1, 0, 1)) / 3;
    grid(0, dim(1) + 1, 0) = (grid(1, dim(1) + 1, 0) + grid(0, dim(1), 0) +
                              grid(0, dim(1) + 1, 1)) / 3;
    grid(0, 0, dim(2) + 1) = (grid(1, 0, dim(2) + 1) + grid(0, 1, dim(2) + 1) +
                              grid(0, 0, dim(2))) / 3;
    grid(dim(0) + 1, dim(1) + 1, 0) = (grid(dim(0), dim(1) + 1, 0) +
                                       grid(dim(0) + 1, dim(1), 0) +
                                       grid(dim(0) + 1, dim(1) + 1, 1)) / 3;
    grid(dim(0) + 1, 0, dim(2) + 1) = (grid(dim(0), 0, dim(2) + 1) +
                                       grid(dim(0) + 1, 1, dim(2) + 1) +
                                       grid(dim(0) + 1, 0, dim(2))) / 3;
    grid(0, dim(1) + 1, dim(2) + 1) = (grid(1, dim(1) + 1, dim(2) + 1) +
                                       grid(0, dim(1), dim(2) + 1) +
                                       grid(0, dim(1) + 1, dim(2))) / 3;
    grid(dim(0) + 1, dim(1) + 1, dim(2) + 1) = (grid(dim(0), dim(1) + 1, dim(2) + 1) +
                                                grid(dim(0) + 1, dim(1), dim(2) + 1) +
                                                grid(dim(0) + 1, dim(1) + 1, dim(2))) / 3;
}
void setContinuityBoundaries(FlowGrid &grid, const Indices &dim) {
    setBoundaries(grid, -1, dim);
}
void setHorizontalNeumannBoundaries(FlowGrid &grid, const Indices &dim) {
    setBoundaries(grid, 0, dim);
}
void setVerticalNeumannBoundaries(FlowGrid &grid, const Indices &dim) {
    setBoundaries(grid, 1, dim);
}
void setDepthNeumannBoundaries(FlowGrid &grid, const Indices &dim) {
    setBoundaries(grid, 2, dim);
}
