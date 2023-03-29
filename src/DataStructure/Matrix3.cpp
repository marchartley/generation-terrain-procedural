#include "DataStructure/Matrix3.h"
//#include "Utils/stb_image.h"

template<>
Matrix3<Vector3> Matrix3<Vector3>::curl() {
    Matrix3<Vector3> returningGrid(this->sizeX, this->sizeY, this->sizeZ);
    #pragma omp parallel for collapse(3)
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
    #pragma omp parallel for collapse(3)
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
Vector3 Matrix3<Vector3>::gradient(Vector3 position)
{
    this->raiseErrorOnBadCoord = false;
    Vector3 flooredPos = position.floor();
    Vector3 offset = position - flooredPos;
    return Vector3(
                at(flooredPos + Vector3(1, 0, 0)).x * (1 - offset.x) + at(flooredPos).x * offset.x,
                at(flooredPos + Vector3(0, 1, 0)).y * (1 - offset.y) + at(flooredPos).y * offset.y,
                at(flooredPos + Vector3(0, 0, 1)).z * (1 - offset.z) + at(flooredPos).z * offset.z
                );
}

template<>
Vector3 Matrix3<Vector3>::gradient(float posX, float posY, float posZ)
{
    return gradient(Vector3(posX, posY, posZ));
}

template<>
Matrix3<Vector3> Matrix3<Vector3>::gradient()
{
    this->raiseErrorOnBadCoord = false;
    Matrix3<Vector3> returningGrid(this->sizeX, this->sizeY, this->sizeZ);
    #pragma omp parallel for collapse(3)
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
//                returningGrid.at(x, y, z) = gradient(x, y, z);
//                continue;
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

template<>
Matrix3<Vector3> Matrix3<Vector3>::fromImageRGB(std::string filename)
{
    int imgW, imgH, nbChannels;
    unsigned char *c_data = stbi_load(filename.c_str(), &imgW, &imgH, &nbChannels, STBI_rgb); // Load image, force 3 channel
    if (c_data == NULL)
    {
        std::cerr << "Error : Impossible to load " << filename << "\n";
        std::cerr << "Either file is not found, or type is incorrect. Available file types are : \n";
        std::cerr << "\t- JPG, \n\t- PNG, \n\t- TGA, \n\t- BMP, \n\t- PSD, \n\t- GIF, \n\t- HDR, \n\t- PIC";
        exit (-1);
        return Matrix3<Vector3>();
    }
    float *data = new float[imgW * imgH * 3];
    for (int i = 0; i < imgW * imgH * 3; i++)
        data[i] = c_data[i];
    stbi_image_free(c_data);

    Matrix3<Vector3> map(imgW, imgH);
    for (int x = 0; x < imgW; x++) {
        for (int y = 0; y < imgH; y++) {
            int index = x + y * imgW;
            map.at(x, y) = Vector3(data[3 * index + 0],
                    data[3 * index + 1],
                    data[3 * index + 1]) / 255.f;
        }
    }
    if (data != nullptr)
        delete[] data;//stbi_image_free(data);

    return map;
}

template<>
Matrix3<float> Matrix3<float>::fromImageBW(std::string filename)
{
    int imgW, imgH, nbChannels;
    unsigned char *c_data = stbi_load(filename.c_str(), &imgW, &imgH, &nbChannels, STBI_grey); // Load image, force 1 channel
    if (c_data == NULL)
    {
        std::cerr << "Error : Impossible to load " << filename << "\n";
        std::cerr << "Either file is not found, or type is incorrect. Available file types are : \n";
        std::cerr << "\t- JPG, \n\t- PNG, \n\t- TGA, \n\t- BMP, \n\t- PSD, \n\t- GIF, \n\t- HDR, \n\t- PIC";
        exit (-1);
        return Matrix3<float>();
    }
    float *data = new float[imgW * imgH];
    for (int i = 0; i < imgW * imgH; i++)
        data[i] = c_data[i];
    stbi_image_free(c_data);

    Matrix3<float> map(imgW, imgH);
    for (int x = 0; x < imgW; x++) {
        for (int y = 0; y < imgH; y++) {
            int index = x + y * imgW;
            map.at(x, y) = data[index] / 255.f;
        }
    }
    if (data != nullptr)
        delete[] data;//stbi_image_free(data);

    return map;
}
