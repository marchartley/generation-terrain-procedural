#include "DataStructure/Matrix3.h"


template<>
Matrix3<Vector3> Matrix3<Vector3>::curl() {
    Matrix3<Vector3> returningGrid(this->sizeX, this->sizeY, this->sizeZ);
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                Vector3& vec = this->at(x, y, z);
                returningGrid.at(x, y, z) = Vector3(vec.z - vec.y, vec.x - vec.z, vec.y - vec.x);
            }
        }
    }
    return returningGrid;
}
template<>
Matrix3<Vector3> Matrix3<Vector3>::rot() { return this->curl(); }

template<>
Matrix3<float> Matrix3<Vector3>::divergence()
{
    this->raiseErrorOnBadCoord = false;
    Matrix3<float> returningGrid(this->sizeX, this->sizeY, this->sizeZ);
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                // Need to change the divergence function...
//                returningGrid.at(x, y, z) = this->at(x, y, z).divergence();
                returningGrid.at(x, y, z) = ((this->at(x + 1, y, z) - this->at(x - 1, y, z)).x +
                                             (this->at(x, y + 1, z) - this->at(x, y - 1, z)).y +
                                             (this->at(x, y, z + 1) - this->at(x, y, z - 1)).z) * .5f;
            }
        }
    }
    this->raiseErrorOnBadCoord = true;
    return returningGrid;
}

template<>
Matrix3<Vector3> Matrix3<Vector3>::gradient()
{
    this->raiseErrorOnBadCoord = false;
    Matrix3<Vector3> returningGrid(this->sizeX, this->sizeY, this->sizeZ);
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                // Need to change the divergence function...
//                returningGrid.at(x, y, z) = this->at(x, y, z).divergence();
                returningGrid.at(x, y, z) = Vector3((this->at(x + 1, y, z) - this->at(x - 1, y, z)).x * .5f,
                                                    (this->at(x, y + 1, z) - this->at(x, y - 1, z)).y * .5f,
                                                    (this->at(x, y, z + 1) - this->at(x, y, z - 1)).z * .5f);
            }
        }
    }
    this->raiseErrorOnBadCoord = true;
    return returningGrid;
}

template<>
Matrix3<Vector3> Matrix3<Vector3>::random(size_t sizeX, size_t sizeY, size_t sizeZ)
{
    Matrix3<Vector3> mat(sizeX, sizeY, sizeZ);
    for (Vector3& val : mat.data)
        val = Vector3::random();
    return mat;
}
