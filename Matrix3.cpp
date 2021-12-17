#include "Matrix3.h"


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
    Matrix3<float> returningGrid(this->sizeX, this->sizeY, this->sizeZ);
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                returningGrid.at(x, y, z) = this->at(x, y, z).divergence();
            }
        }
    }
    return returningGrid;
}

template<>
Matrix3<Vector3> Matrix3<Vector3>::random(int sizeX, int sizeY, int sizeZ)
{
    Matrix3<Vector3> mat(sizeX, sizeY, sizeZ);
    for (Vector3& val : mat.data)
        val = Vector3::random();
    return mat;
}
