#ifndef MATH_H
#define MATH_H

#include <array>

#include <unsupported/Eigen/CXX11/Tensor>

typedef float Scalar;
const Eigen::Index kGridDimensions = 3;
typedef Eigen::Tensor<Scalar, kGridDimensions> FlowGrid;
typedef Eigen::Array<Scalar, kGridDimensions, 1> Location;
typedef Eigen::Array<FlowGrid::Index, kGridDimensions, 1> Indices;
typedef std::array<FlowGrid::Index, kGridDimensions> TensorIndices;
typedef std::function<void(FlowGrid&)> BoundarySetter;

void linearSolve(FlowGrid &solution, const FlowGrid &initial, Scalar alpha, Scalar beta,
                 const Indices &dim, BoundarySetter setBoundaries,
                 unsigned int iterations = 20);

// Linearly interpolates grid to nearest neighbors
Scalar interpolate(const FlowGrid &grid, Location x);

// Sets boundary conditions on grids
void setBoundaries(FlowGrid &grid, int b, const Indices &dim);
void setContinuityBoundaries(FlowGrid &grid, const Indices &dim);
void setHorizontalNeumannBoundaries(FlowGrid &grid, const Indices &dim);
void setVerticalNeumannBoundaries(FlowGrid &grid, const Indices &dim);
void setDepthNeumannBoundaries(FlowGrid &grid, const Indices &dim);

#endif // MATH_H
