#include "Utils/Globals.h"
#include "TerrainGen/Grid.h"

#include "Utils/FastNoiseLit.h"
#include "Utils/Utils.h"

#define STB_IMAGE_IMPLEMENTATION
#include "Utils/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "Utils/stb_image_write.h"


Grid::Grid(int nx, int ny, float maxHeight, float tileSize)
    : maxHeight(maxHeight), tileSize(tileSize) {
    // Create and configure FastNoise object
    FastNoiseLite noise;
    noise.SetFrequency(0.01);
    noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
    noise.SetFractalType(FastNoiseLite::FractalType_FBm);
    noise.SetFractalLacunarity(2.0);
    noise.SetFractalGain(0.7);
    noise.SetFractalWeightedStrength(0.5);
    noise.SetFractalOctaves(10);

    // Make a first pass, to get some noise
    float min = 1000, max = -1000;
    this->vertices = Matrix3<Vector3>(nx, ny);
    this->normals = Matrix3<Vector3>(nx, ny);
    for (int x = 0; x < this->getSizeX(); x++) {
        for (int y = 0; y < this->getSizeY(); y++) {
            float z = noise.GetNoise((float)x, (float)y);
            this->vertices.at(x, y) = Vector3(x, y, z);
            min = min < z ? min : z;
            max = max > z ? max : z;
        }
    }
    // Sadely, we need a second pass because the noise function is a little f*cked up... we don't always have values in [-1;1)
    for (int x = 0; x < this->getSizeX(); x++) {
        for (int y = 0; y < this->getSizeY(); y++) {
            this->vertices.at(x, y).z = maxHeight * (this->vertices.at(x, y).z - min) / (max - min);
        }
    }
    this->computeNormals();
}

Grid::Grid(std::string heightmap_filename, int nx, int ny, float max_height, float tileSize)
    : maxHeight(max_height), tileSize(tileSize)
{
    this->loadFromHeightmap(heightmap_filename, nx, ny, max_height, tileSize);
}

Grid::Grid() : Grid(10, 10, 5.0) {

}

void Grid::createMesh()
{
    std::vector<Vector3> vecs;

    for(int x = 0; x < this->getSizeX() - 1; x++) {
        for (int y = 0; y < this->getSizeY() - 1; y++) {
            vecs.push_back(this->vertices.at(x, y));
            vecs.push_back(this->vertices.at(x, y+1));
            vecs.push_back(this->vertices.at(x+1, y+1));

            vecs.push_back(this->vertices.at(x, y));
            vecs.push_back(this->vertices.at(x+1, y+1));
            vecs.push_back(this->vertices.at(x+1, y));
        }
    }
    this->mesh.fromArray(vecs);
//    return this->vertexArrayFloat;
}

void Grid::computeNormals() {
//    Vector3 upVector(0.0, 0.0, 1.0);
//    for (int x = 0; x < this->sizeX; x++)
//        for (int y = 0; y < this->sizeY; y++)
//            this->normals.at(x, y) = upVector;

    this->normals = Matrix3<Vector3>(this->getSizeX(), this->getSizeY());
    this->vertices.raiseErrorOnBadCoord = false;
    for (int x = 0; x < this->getSizeX(); x++) {
        for (int y = 0; y < this->getSizeY(); y++) {
            Vector3 pos = this->vertices.at(x, y);
            Vector3 v1 = this->vertices.at(x - 1, y    ) - pos;
            Vector3 v2 = this->vertices.at(x - 1, y - 1) - pos;
            Vector3 v3 = this->vertices.at(x    , y - 1) - pos;
            Vector3 v4 = this->vertices.at(x + 1, y - 1) - pos;
            Vector3 v5 = this->vertices.at(x + 1, y    ) - pos;
            Vector3 v6 = this->vertices.at(x + 1, y + 1) - pos;
            Vector3 v7 = this->vertices.at(x    , y + 1) - pos;
            Vector3 v8 = this->vertices.at(x - 1, y + 1) - pos;

            this->normals.at(x, y) += v1.cross(v2);
            this->normals.at(x, y) += v2.cross(v3);
            this->normals.at(x, y) += v3.cross(v4);
            this->normals.at(x, y) += v4.cross(v5);
            this->normals.at(x, y) += v5.cross(v6);
            this->normals.at(x, y) += v6.cross(v7);
            this->normals.at(x, y) += v7.cross(v8);
            this->normals.at(x, y) += v8.cross(v1);
            this->normals.at(x, y).normalize();
        }
    }
}

void Grid::display(bool displayNormals) {
    this->mesh.display();
    /*
    glPushMatrix();

    float maxi = -9999999, mini = 9999999;
    for (int x = 0; x < this->sizeX; x++)
        for(int y = 0; y < this->sizeY; y++){
            if (this->vertices.at(x, y).z > maxi)
                maxi = this->vertices.at(x, y).z;
            if (this->vertices.at(x, y).z < mini)
                mini = this->vertices.at(x, y).z;
        }
    glScalef(1/this->tileSize, 1/this->tileSize, 1/this->tileSize);
    glTranslatef(-(sizeX - 1) / 2.0, - (sizeY - 1) / 2.0, - (maxi + mini) / 2.0);

    glBegin(GL_TRIANGLES);
//    glColor4f(1.0, 0.5, 0.5, 0.0);
    for(int x = 0; x < this->sizeX - 1; x++) {
        for (int y = 0; y < this->sizeY - 1; y++) {
            glVertex3f(this->vertices.at(x, y).x, this->vertices.at(x, y).y, this->vertices.at(x, y).z);
            glVertex3f(this->vertices.at(x, y+1).x, this->vertices.at(x, y+1).y, this->vertices.at(x, y+1).z);
            glVertex3f(this->vertices.at(x+1, y+1).x, this->vertices.at(x+1, y+1).y, this->vertices.at(x+1, y+1).z);

            glVertex3f(this->vertices.at(x, y).x, this->vertices.at(x, y).y, this->vertices.at(x, y).z);
            glVertex3f(this->vertices.at(x+1, y).x, this->vertices.at(x+1, y).y, this->vertices.at(x+1, y).z);
            glVertex3f(this->vertices.at(x+1, y+1).x, this->vertices.at(x+1, y+1).y, this->vertices.at(x+1, y+1).z);
        }
    }
    glEnd();

    if (displayNormals) {
        this->computeNormals();
        for(int x = 0; x < this->sizeX - 1; x++) {
            for (int y = 0; y < this->sizeY - 1; y++) {
                glBegin(GL_LINES);
                glVertex3f(this->vertices.at(x, y).x, this->vertices.at(x, y).y, this->vertices.at(x, y).z);
                glVertex3f(this->vertices.at(x, y).x + this->normals.at(x, y).x, this->vertices.at(x, y).y + this->normals.at(x, y).y, this->vertices.at(x, y).z + this->normals.at(x, y).z);
                glEnd();
            }
        }
    }
    glPopMatrix();*/
}

float Grid::getMaxHeight()
{
    float max = 0;
    for (const Vector3& v : this->vertices)
        max = std::max(max, v.z);
    return max;
}

void Grid::fromVoxelGrid(VoxelGrid &voxelGrid) {
    this->vertices = Matrix3<Vector3>(voxelGrid.getSizeX(), voxelGrid.getSizeY());
    for (int x = 0; x < voxelGrid.getSizeX(); x++) {
        for (int y = 0; y < voxelGrid.getSizeY(); y++) {
            this->vertices.at(x, y) = Vector3(x, y, voxelGrid.getHeight(x, y)/*/ (float)(voxelGrid.getSizeZ())*/);
        }
    }
    this->computeNormals();
    this->createMesh();
}

void Grid::loadFromHeightmap(std::string heightmap_filename, int nx, int ny, float max_height, float tileSize)
{
    this->tileSize = tileSize;
    this->maxHeight = max_height;
    int imgW, imgH, nbChannels;
    unsigned char *data = stbi_load(heightmap_filename.c_str(), &imgW, &imgH, &nbChannels, STBI_grey); // Load image, force 1 channel
    if (data == NULL)
    {
        std::cerr << "Error : Impossible to load " << heightmap_filename << "\n";
        std::cerr << "Either file is not found, or type is incorrect. Available file types are : \n";
        std::cerr << "\t- JPG, \n\t- PNG, \n\t- TGA, \n\t- BMP, \n\t- PSD, \n\t- GIF, \n\t- HDR, \n\t- PIC";
        exit (-1);
        return;
    }
    if (nx == -1)
        nx = imgW;
    if (ny == -1)
        ny = imgH;

    float max = 0;

    Matrix3<float> map(imgW, imgH);
    for (int x = 0; x < imgW; x++) {
        for (int y = 0; y < imgH; y++) {
            float value = (float)data[x + y * imgW];
            map.at(x, y) = value;
            max = std::max(max, value);
        }
    }
    stbi_image_free(data);

    map = map.resize(nx, ny, 1);
    if (this->maxHeight == -1) {
        maxHeight = max;
    } else {
        map *= (maxHeight / max);
    }
    this->vertices = Matrix3<Vector3>(map.sizeX, map.sizeY);
    this->normals = Matrix3<Vector3>(map.sizeX, map.sizeY);
    for (int x = 0; x < map.sizeX; x++) {
        for (int y = 0; y < map.sizeY; y++) {
            this->vertices.at(x, y) = Vector3(x, y, map.at(x, y));
        }
    }
    this->computeNormals();
}

void Grid::saveHeightmap(std::string heightmap_filename)
{
    std::string ext = toUpper(getExtention(heightmap_filename));
    int width = this->getSizeX();
    int height = this->getSizeY();
    // To heightmap
    std::vector<float> toFloatData(width*height);
    std::vector<uint8_t> toIntData(width*height);

    float newHeight = this->maxHeight;
    for (const auto& vert : this->vertices)
        newHeight = std::max(newHeight, vert.z);

    for (size_t i = 0; i < this->vertices.size(); i++) {
        toFloatData[i] = this->vertices[i].z / newHeight;
        toIntData[i] = toFloatData[i] * 255;
    }
    if (ext == "PNG")
        stbi_write_png(heightmap_filename.c_str(), this->getSizeX(), this->getSizeY(), 1, toIntData.data(), this->getSizeX() * 1);
    else if (ext == "JPG")
        stbi_write_jpg(heightmap_filename.c_str(), this->getSizeX(), this->getSizeY(), 1, toIntData.data(), 95);
    else if (ext == "BMP")
        stbi_write_bmp(heightmap_filename.c_str(), this->getSizeX(), this->getSizeY(), 1, toIntData.data());
    else if (ext == "TGA")
        stbi_write_tga(heightmap_filename.c_str(), this->getSizeX(), this->getSizeY(), 1, toIntData.data());
    else if (ext == "HDR")
        stbi_write_hdr(heightmap_filename.c_str(), this->getSizeX(), this->getSizeY(), 1, toFloatData.data());
    else {
        std::cerr << "Trying to save map without valid extension. Possible extensions :\n\t- png\n\t- jpg\n\t- tga\n\t- bmp\n\t- hdr" << std::endl;
    }
}
