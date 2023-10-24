#ifndef VECTORFIELD_H
#define VECTORFIELD_H

#include "math.h"

template<FlowGrid::Index numStaggers, std::size_t numCoords>
class VectorField {
public:
    VectorField(const TensorIndices &dimensions);

    static const std::size_t coords = numCoords;

    void clear();

    const FlowGrid &operator[](std::size_t coord) const;
    FlowGrid &operator[](std::size_t coord);

    VectorField<numStaggers, numCoords>
    &operator+=(const VectorField<numStaggers, numCoords> &rhs);
    VectorField<numStaggers, numCoords>
    &operator-=(const VectorField<numStaggers, numCoords> &rhs);
    VectorField<numStaggers, numCoords> &operator*=(Scalar rhs);

private:
    std::array<FlowGrid, numCoords> grids;

};
template<FlowGrid::Index numStaggers, std::size_t numCoords>
VectorField<numStaggers, numCoords>
operator-=(VectorField<numStaggers, numCoords> lhs,
           const VectorField<numStaggers, numCoords> &rhs);
template<FlowGrid::Index numStaggers, std::size_t numCoords>
VectorField<numStaggers, numCoords>
operator-(VectorField<numStaggers, numCoords> lhs,
          const VectorField<numStaggers, numCoords> &rhs);
template<FlowGrid::Index numStaggers, std::size_t numCoords>
VectorField<numStaggers, numCoords>
operator*(VectorField<numStaggers, numCoords> lhs, Scalar rhs);
template<FlowGrid::Index numStaggers, std::size_t numCoords>
VectorField<numStaggers, numCoords>
operator*(Scalar lhs, VectorField<numStaggers, numCoords> rhs);

#include "vectorfield.tpp"

#endif // VECTORFIELD_H
