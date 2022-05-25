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
//    this->vertices = Matrix3<Vector3>(nx, ny);
    this->normals = Matrix3<Vector3>(nx, ny);
    this->heights = Matrix3<float>(nx, ny);
    for (int x = 0; x < this->getSizeX(); x++) {
        for (int y = 0; y < this->getSizeY(); y++) {
            float z = noise.GetNoise((float)x, (float)y);
            this->heights.at(x, y) = z;
//            this->vertices.at(x, y) = Vector3(x, y, z);
            min = min < z ? min : z;
            max = max > z ? max : z;
        }
    }
    this->heights = ((this->heights - min) / (max - min)) * maxHeight;
    // Sadely, we need a second pass because the noise function is a little f*cked up... we don't always have values in [-1;1)
//    for (int x = 0; x < this->getSizeX(); x++) {
//        for (int y = 0; y < this->getSizeY(); y++) {
//            this->vertices.at(x, y).z = maxHeight * (this->vertices.at(x, y).z - min) / (max - min);
//        }
//    }
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
            vecs.push_back(Vector3(x  , y  , this->heights.at(x  , y  )));
            vecs.push_back(Vector3(x  , y+1, this->heights.at(x  , y+1)));
            vecs.push_back(Vector3(x+1, y+1, this->heights.at(x+1, y+1)));

            vecs.push_back(Vector3(x  , y  , this->heights.at(x  , y  )));
            vecs.push_back(Vector3(x+1, y+1, this->heights.at(x+1, y+1)));
            vecs.push_back(Vector3(x+1, y  , this->heights.at(x+1, y  )));
//            vecs.push_back(this->vertices.at(x, y));
//            vecs.push_back(this->vertices.at(x, y+1));
//            vecs.push_back(this->vertices.at(x+1, y+1));

//            vecs.push_back(this->vertices.at(x, y));
//            vecs.push_back(this->vertices.at(x+1, y+1));
//            vecs.push_back(this->vertices.at(x+1, y));
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
    this->heights.raiseErrorOnBadCoord = false;
//    this->vertices.raiseErrorOnBadCoord = false;
    for (int x = 0; x < this->getSizeX(); x++) {
        for (int y = 0; y < this->getSizeY(); y++) {
            /*Vector3 pos = this->vertices.at(x, y);
            Vector3 v1 = this->vertices.at(x - 1, y    ) - pos;
            Vector3 v2 = this->vertices.at(x - 1, y - 1) - pos;
            Vector3 v3 = this->vertices.at(x    , y - 1) - pos;
            Vector3 v4 = this->vertices.at(x + 1, y - 1) - pos;
            Vector3 v5 = this->vertices.at(x + 1, y    ) - pos;
            Vector3 v6 = this->vertices.at(x + 1, y + 1) - pos;
            Vector3 v7 = this->vertices.at(x    , y + 1) - pos;
            Vector3 v8 = this->vertices.at(x - 1, y + 1) - pos;
            */
            Vector3 pos = Vector3(x   , y    , this->heights.at(x, y));
            Vector3 v1 = Vector3(x - 1, y    , this->heights.at(x - 1, y    )) - pos;
            Vector3 v2 = Vector3(x - 1, y - 1, this->heights.at(x - 1, y - 1)) - pos;
            Vector3 v3 = Vector3(x    , y - 1, this->heights.at(x    , y - 1)) - pos;
            Vector3 v4 = Vector3(x + 1, y - 1, this->heights.at(x + 1, y - 1)) - pos;
            Vector3 v5 = Vector3(x + 1, y    , this->heights.at(x + 1, y    )) - pos;
            Vector3 v6 = Vector3(x + 1, y + 1, this->heights.at(x + 1, y + 1)) - pos;
            Vector3 v7 = Vector3(x    , y + 1, this->heights.at(x    , y + 1)) - pos;
            Vector3 v8 = Vector3(x - 1, y + 1, this->heights.at(x - 1, y + 1)) - pos;

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
    return this->heights.max();
//    float max = 0;
//    for (const Vector3& v : this->vertices)
//        max = std::max(max, v.z);
//    return max;
}

std::vector<std::vector<Vector3>> Grid::hydraulicErosion()
{
    std::vector<std::vector<Vector3>> traces;
    int nb_particles = 1000;

//    int erosionRadius = 1;
    float inertia = 0; //.05f;
//    float sedimentCapacityFactor = 4;
    float minSedimentCapacity = .01f;
    float minWater = .01f;
//    float maxCapacity = 1.f;
    float depositRate = .5f;
    float erosionRate = .5f;
//    float depositRate = .3f;
//    float erosionRate = .3f;
    float evaporationRate = .001f;
    float initialWaterVolume = 1;
    float initialSpeed = 1;
//    float gravity = 4;
    float dt = 1.2;

    float friction = .05f;

//    Matrix3<float> heights(vertices.getDimensions());
//    for (size_t i = 0; i < vertices.size(); i++)
//        heights[i] = vertices[i].z;
//    heights /= this->getMaxHeight();
    Matrix3<Vector3> gradient;

    int maxIterationPerParticle = 300;

    for (int i = 0; i < nb_particles; i++) {


        std::vector<Vector3> trace;
        gradient = heights.gradient() * -1.f;
        for (int x = 0; x < gradient.sizeX; x++) {
            if (x == 0 || x == gradient.sizeX - 1)
                for (int y = 0; y < gradient.sizeY; y++)
                    gradient.at(x, y) = Vector3();
            gradient.at(x, 0) = Vector3();
            gradient.at(x, gradient.sizeY - 1) = Vector3();
        }
        Vector3 position = Vector3(random_gen::generate(1, this->getSizeX() - 1), random_gen::generate(1, this->getSizeY() - 1));

        Vector3 direction = gradient.interpolate(position) * dt; //.normalize();
        float speed = initialSpeed;
        float water = initialWaterVolume;
        float sediments = 0;

        Vector3 meanPos = position;
        int iteration = 0;
        while (heights.checkCoord(position) && iteration < maxIterationPerParticle && water > minWater) {
            float height = heights.interpolate(position);

            direction += (gradient.interpolate(position) / water) * dt;
            position += direction * dt;
            direction *= (1.0 - dt * friction);

            float newHeight = heights.interpolate(position);
            float deltaHeight = newHeight - height;

            trace.push_back(Vector3(position.x, position.y, heights.interpolate(position.x, position.y, 0) * this->getMaxHeight()));
            if (!heights.checkCoord(position + direction * 3.f)) break;

            float maxSediments = std::max(0.f, water * direction.norm() * deltaHeight);
            float diff = maxSediments - sediments;
            sediments += dt * depositRate * diff;

            float erosion = dt * depositRate * diff * water;

            float xOffset = position.x - position.rounded().x;
            float yOffset = position.y - position.rounded().y;
            if (heights.checkCoord(position.rounded() + Vector3(0, 0, 0)))
                heights.at(position.rounded() + Vector3(0, 0, 0)) -= erosion * (1 - xOffset) * (1 - yOffset);
            if (heights.checkCoord(position.rounded() + Vector3(1, 0, 0)))
                heights.at(position.rounded() + Vector3(1, 0, 0)) -= erosion * (xOffset) * (1 - yOffset);
            if (heights.checkCoord(position.rounded() + Vector3(0, 1, 0)))
                heights.at(position.rounded() + Vector3(0, 1, 0)) -= erosion * (1 - xOffset) * (yOffset);
            if (heights.checkCoord(position.rounded() + Vector3(1, 1, 0)))
                heights.at(position.rounded() + Vector3(1, 1, 0)) -= erosion * (xOffset) * (yOffset);

            water *= (1.0 - dt * evaporationRate);
/*
            float height = heights.interpolate(position);
            direction = (direction * inertia + (gradient.interpolate(position) * (1.f - inertia))).normalize();
            position += direction;
//            meanPos += position;
//            if ((meanPos / (float)iteration - position).norm2() < 0.5f) break;
            trace.push_back(Vector3(position.x, position.y, heights.interpolate(position.x, position.y, 0) * this->getMaxHeight()));
            if (!heights.checkCoord(position + direction * 3.f)) break;

            float newHeight = heights.interpolate(position);
            float deltaHeight = newHeight - height;

            float erosion = erosionRate * gradient.interpolate(position).norm(); // * deltaHeight;// * erosionRate;
            float deposit = sediments * depositRate;

            deposit -= erosion;
            sediments -= deposit;
            float xOffset = position.x - position.rounded().x;
            float yOffset = position.y - position.rounded().y;
            if (heights.checkCoord(position.rounded() + Vector3(0, 0, 0)))
                heights.at(position.rounded() + Vector3(0, 0, 0)) += deposit * (1 - xOffset) * (1 - yOffset);
            if (heights.checkCoord(position.rounded() + Vector3(1, 0, 0)))
                heights.at(position.rounded() + Vector3(1, 0, 0)) += deposit * (xOffset) * (1 - yOffset);
            if (heights.checkCoord(position.rounded() + Vector3(0, 1, 0)))
                heights.at(position.rounded() + Vector3(0, 1, 0)) += deposit * (1 - xOffset) * (yOffset);
            if (heights.checkCoord(position.rounded() + Vector3(1, 1, 0)))
                heights.at(position.rounded() + Vector3(1, 1, 0)) += deposit * (xOffset) * (yOffset);
*/
            /*

            direction = (direction * inertia + (gradient.interpolate(position) * (1.f - inertia))).normalize();
            position += direction;
            if (!heights.checkCoord(position)) break;

            float newHeight = heights.interpolate(position);
            float deltaHeight = newHeight - height;
            trace.push_back(Vector3(position.x, position.y, getHeight(position.x, position.y)));
//            float speed = direction.norm();
            float sedimentsCapacity = std::max(-deltaHeight * speed * water * sedimentCapacityFactor, minSedimentCapacity);

            if (sediments > sedimentsCapacity || deltaHeight > 0) {
                // Deposit
                float deposit = (deltaHeight > 0) ? std::min(deltaHeight, sediments) : (sediments - sedimentsCapacity) * depositRate;
                sediments -= deposit;
                if (sediments < -0.5)
                    std::cout << "Deposit way too big..." << std::endl;

                float xOffset = position.x - position.rounded().x;
                float yOffset = position.y - position.rounded().y;
                if (heights.checkCoord(position.rounded() + Vector3(0, 0, 0)))
                    heights.at(position.rounded() + Vector3(0, 0, 0)) += deposit * (1 - xOffset) * (1 - yOffset);
                if (heights.checkCoord(position.rounded() + Vector3(1, 0, 0)))
                    heights.at(position.rounded() + Vector3(1, 0, 0)) += deposit * (xOffset) * (1 - yOffset);
                if (heights.checkCoord(position.rounded() + Vector3(0, 1, 0)))
                    heights.at(position.rounded() + Vector3(0, 1, 0)) += deposit * (1 - xOffset) * (yOffset);
                if (heights.checkCoord(position.rounded() + Vector3(1, 1, 0)))
                    heights.at(position.rounded() + Vector3(1, 1, 0)) += deposit * (xOffset) * (yOffset);
            } else {
                // Erosion
                float erosion = std::min(heights.interpolate(position), std::min((sedimentsCapacity - sediments) * erosionRate, -deltaHeight));
                sediments += erosion;

                float xOffset = position.x - position.rounded().x;
                float yOffset = position.y - position.rounded().y;
                if (heights.checkCoord(position.rounded() + Vector3(0, 0, 0)))
                    heights.at(position.rounded() + Vector3(0, 0, 0)) -= erosion * (1 - xOffset) * (1 - yOffset);
                if (heights.checkCoord(position.rounded() + Vector3(1, 0, 0)))
                    heights.at(position.rounded() + Vector3(1, 0, 0)) -= erosion * (xOffset) * (1 - yOffset);
                if (heights.checkCoord(position.rounded() + Vector3(0, 1, 0)))
                    heights.at(position.rounded() + Vector3(0, 1, 0)) -= erosion * (1 - xOffset) * (yOffset);
                if (heights.checkCoord(position.rounded() + Vector3(1, 1, 0)))
                    heights.at(position.rounded() + Vector3(1, 1, 0)) -= erosion * (xOffset) * (yOffset);
            }
            speed = std::sqrt(speed * speed + deltaHeight * gravity);
            water *= (1 - evaporationRate);
            */
            iteration ++;
        }
        traces.push_back(trace);
    }
//    heights *= this->getMaxHeight();
//    for (size_t i = 0; i < vertices.size(); i++)
//        vertices[i].z = heights[i] ;

    this->computeNormals();
    this->createMesh();

    return traces;
}

void Grid::thermalErosion()
{

}

void Grid::windErosion()
{

}

void Grid::fromVoxelGrid(VoxelGrid &voxelGrid) {
    this->vertices = Matrix3<Vector3>(voxelGrid.getSizeX(), voxelGrid.getSizeY());
    for (int x = 0; x < voxelGrid.getSizeX(); x++) {
        for (int y = 0; y < voxelGrid.getSizeY(); y++) {
            this->heights.at(x, y) = voxelGrid.getHeight(x, y);
//            this->vertices.at(x, y) = Vector3(x, y, voxelGrid.getHeight(x, y)/*/ (float)(voxelGrid.getSizeZ())*/);
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
    this->heights = map;
    /*
    this->vertices = Matrix3<Vector3>(map.sizeX, map.sizeY);
    this->normals = Matrix3<Vector3>(map.sizeX, map.sizeY);
    for (int x = 0; x < map.sizeX; x++) {
        for (int y = 0; y < map.sizeY; y++) {
            this->vertices.at(x, y) = Vector3(x, y, map.at(x, y));
        }
    }*/
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

    float newHeight = std::max(this->maxHeight, this->heights.max());
    /*
    float newHeight = this->maxHeight;
    for (const auto& vert : this->vertices)
        newHeight = std::max(newHeight, vert.z);
    */
    toFloatData = (this->heights/newHeight).data;
    for (size_t i = 0; i < this->heights.size(); i++) {
//        toFloatData[i] = this->vertices[i].z / newHeight;
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
