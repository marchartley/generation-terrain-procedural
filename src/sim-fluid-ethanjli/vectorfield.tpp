#include "vectorfield.h"

template<FlowGrid::Index numStaggers, std::size_t numCoords>
VectorField<numStaggers, numCoords>::VectorField(const TensorIndices &dimensions) {
    for (auto &grid : grids) {
      grid = FlowGrid(dimensions);
    }
    clear();
}
template<FlowGrid::Index numStaggers, std::size_t numCoords>
void VectorField<numStaggers, numCoords>::clear() {
    for (auto &grid : grids) {
        grid.setZero();
    }
}
template<FlowGrid::Index numStaggers, std::size_t numCoords>
const FlowGrid &VectorField<numStaggers, numCoords>::operator[](std::size_t coord) const {
    return grids[coord];
}
template<FlowGrid::Index numStaggers, std::size_t numCoords>
FlowGrid &VectorField<numStaggers, numCoords>::operator[](std::size_t coord) {
    return grids[coord];
}
template<FlowGrid::Index numStaggers, std::size_t numCoords>
VectorField<numStaggers, numCoords>
&VectorField<numStaggers, numCoords>::operator+=(const VectorField<numStaggers, numCoords> &rhs) {
    for (std::size_t i = 0; i < numCoords; ++i) {
        grids[i] += rhs.grids[i];
    }
    return *this;
}
template<FlowGrid::Index numStaggers, std::size_t numCoords>
VectorField<numStaggers, numCoords>
&VectorField<numStaggers, numCoords>::operator-=(const VectorField<numStaggers, numCoords> &rhs) {
    for (std::size_t i = 0; i < numCoords; ++i) {
        grids[i] -= rhs.grids[i];
    }
    return *this;
}
template<FlowGrid::Index numStaggers, std::size_t numCoords>
VectorField<numStaggers, numCoords>
&VectorField<numStaggers, numCoords>::operator*=(Scalar rhs) {
    for (std::size_t i = 0; i < numCoords; ++i) {
        grids[i] = grids[i] * rhs;
    }
    return *this;
}

template<FlowGrid::Index numStaggers, std::size_t numCoords>
VectorField<numStaggers, numCoords>
operator+(VectorField<numStaggers, numCoords> lhs,
          const VectorField<numStaggers, numCoords> &rhs) {
    lhs += rhs;
    return lhs;
}
template<FlowGrid::Index numStaggers, std::size_t numCoords>
VectorField<numStaggers, numCoords>
operator-(VectorField<numStaggers, numCoords> lhs,
          const VectorField<numStaggers, numCoords> &rhs) {
    lhs -= rhs;
    return lhs;
}
template<FlowGrid::Index numStaggers, std::size_t numCoords>
VectorField<numStaggers, numCoords>
operator*(VectorField<numStaggers, numCoords> lhs, Scalar rhs) {
    lhs *= rhs;
    return lhs;
}
template<FlowGrid::Index numStaggers, std::size_t numCoords>
VectorField<numStaggers, numCoords>
operator*(Scalar lhs, VectorField<numStaggers, numCoords> rhs) {
    rhs *= lhs;
    return rhs;
}
