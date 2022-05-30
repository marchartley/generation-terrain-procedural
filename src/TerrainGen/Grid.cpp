#include "Utils/Globals.h"
#include "TerrainGen/Grid.h"

#include "Utils/FastNoiseLit.h"
#include "Utils/Utils.h"
#include "Utils/BSpline.h"

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
    this->normals = Matrix3<Vector3>(nx, ny);
    this->heights = Matrix3<float>(nx, ny);
    for (int x = 0; x < this->getSizeX(); x++) {
        for (int y = 0; y < this->getSizeY(); y++) {
            float z = noise.GetNoise((float)x, (float)y);
            this->heights.at(x, y) = z;
            min = min < z ? min : z;
            max = max > z ? max : z;
        }
    }
    this->heights = ((this->heights - min) / (max - min)) * maxHeight;
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
    /*
    std::vector<Vector3> vecs;

    for(int x = 0; x < this->getSizeX() - 1; x++) {
        for (int y = 0; y < this->getSizeY() - 1; y++) {
            vecs.push_back(Vector3(x  , y  , this->heights.at(x  , y  )));
            vecs.push_back(Vector3(x  , y+1, this->heights.at(x  , y+1)));
            vecs.push_back(Vector3(x+1, y+1, this->heights.at(x+1, y+1)));

            vecs.push_back(Vector3(x  , y  , this->heights.at(x  , y  )));
            vecs.push_back(Vector3(x+1, y+1, this->heights.at(x+1, y+1)));
            vecs.push_back(Vector3(x+1, y  , this->heights.at(x+1, y  )));
        }
    }
    this->mesh.fromArray(vecs);
    */
}

void Grid::computeNormals() {
    return;
    /*
    this->normals = Matrix3<Vector3>(this->getSizeX(), this->getSizeY());
    this->heights.raiseErrorOnBadCoord = false;

    for (int x = 0; x < this->getSizeX(); x++) {
        for (int y = 0; y < this->getSizeY(); y++) {
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
    }*/
}

void Grid::display(bool displayNormals) {
    this->mesh.display();
}

float Grid::getMaxHeight()
{
    return this->heights.max();
}


// Greatly inspired by Sebastian Lague https://github.com/SebLague/Hydraulic-Erosion
std::vector<std::vector<Vector3>> Grid::hydraulicErosion(bool applyDeposit)
{
    std::vector<std::vector<Vector3>> traces;
    float currentMaxHeight = this->heights.max();
    this->heights /= currentMaxHeight;

    int erosionRadius = 10;
    float inertia = .05f; // At zero, water will instantly change direction to flow downhill. At 1, water will never change direction.
    float sedimentCapacityFactor = 1; // Multiplier for how much sediment a droplet can carry
    float minSedimentCapacity = .01f; // Used to prevent carry capacity getting too close to zero on flatter terrain
    float erodeSpeed = .3f;
    float depositSpeed = .3f;
    float evaporateSpeed = .01f;
    float gravity = 4;
    int maxDropletLifetime = 30;

    float initialWaterVolume = 1;
    float initialSpeed = 1;

    int numIterations = 1000;

    for (int iteration = 0; iteration < numIterations; iteration++) {
        std::vector<Vector3> trace;
        // Create water droplet at random point on map
        Vector3 pos(random_gen::generate(0, this->getSizeX() - 1), random_gen::generate(0, this->getSizeY() - 1), 0);
        Vector3 dir(0, 0, 0);
        float speed = initialSpeed;
        float water = initialWaterVolume;
        float sediment = 0;

        Matrix3<Vector3> gradients = heights.gradient();
        Matrix3<float> precomputedHeights = heights;

        for (int lifetime = 0; lifetime < maxDropletLifetime; lifetime++) {
            trace.push_back(pos);
            if (!applyDeposit)
                sediment = 0;

            // Calculate droplet's height and direction of flow with bilinear interpolation of surrounding heights
            float height = heights.interpolate(pos);
            Vector3 gradient = heights.gradient(pos);

            // Update the droplet's direction and position (move position 1 unit regardless of speed)
            dir = (dir * inertia - gradient * (1 - inertia)).normalize();
            pos += dir;

            // Stop simulating droplet if it's not moving or has flowed over edge of map
            if ((dir.x == 0 && dir.y == 0) || !heights.checkCoord(pos)) {
                break;
            }

            // Find the droplet's new height and calculate the deltaHeight
            float newHeight = heights.interpolate(pos);
            float deltaHeight = newHeight - height;

            // Calculate the droplet's sediment capacity (higher when moving fast down a slope and contains lots of water)
            float sedimentCapacity = std::max (-deltaHeight * speed * water * sedimentCapacityFactor, minSedimentCapacity);

            // If carrying more sediment than capacity, or if flowing uphill:
            if (sediment > sedimentCapacity || deltaHeight > 0) {
                // If moving uphill (deltaHeight > 0) try fill up to the current height, otherwise deposit a fraction of the excess sediment
                float amountToDeposit = (deltaHeight > 0) ? std::min (deltaHeight, sediment) : (sediment - sedimentCapacity) * depositSpeed;
                sediment -= amountToDeposit;

                // Add the sediment to the four nodes of the current cell using bilinear interpolation
                // Deposition is not distributed over a radius (like erosion) so that it can fill small pits
                precomputedHeights.raiseErrorOnBadCoord = false;
                Matrix3<float> brush = Matrix3<float>::gaussian(erosionRadius * 2, erosionRadius * 2, 1, 2.f, (pos - pos.floor())) * amountToDeposit * 1.f;
                for (int x = 0; x < brush.sizeX; x++) {
                    for (int y = 0; y < brush.sizeY; y++) {
                        Vector3 offset(x - erosionRadius, y - erosionRadius);
                        precomputedHeights.raiseErrorOnBadCoord = false;
                        precomputedHeights.at(pos + offset) += brush.at(x, y);
                    }
                }

            } else {
                // Erode a fraction of the droplet's current carry capacity.
                // Clamp the erosion to the change in height so that it doesn't dig a hole in the terrain behind the droplet
                float amountToErode = std::min ((sedimentCapacity - sediment) * erodeSpeed, -deltaHeight);
                Matrix3<float> brush = Matrix3<float>::gaussian(erosionRadius * 2, erosionRadius * 2, 1, 2.f, (pos - pos.floor())) * amountToErode * 1.f;
                // Use erosion brush to erode from all nodes inside the droplet's erosion radius
                for (int x = 0; x < brush.sizeX; x++) {
                    for (int y = 0; y < brush.sizeY; y++) {
                        Vector3 offset(x - erosionRadius, y - erosionRadius);
                        float deltaSediment = std::max(0.f, std::min(heights.interpolate(pos + offset), brush.at(x, y)));
                        precomputedHeights.raiseErrorOnBadCoord = false;
                        precomputedHeights.at(pos + offset) -= deltaSediment;
                        sediment += deltaSediment;
                    }
                }
            }

            // Update droplet's speed and water content
            speed = std::sqrt (std::abs(speed * speed + deltaHeight * gravity));
            water *= (1 - evaporateSpeed);
        }
        heights = precomputedHeights;
        for(auto& h : heights)
            h = std::min(std::max(h, 0.f), 1.f);
        traces.push_back(trace);
    }
    this->heights *= currentMaxHeight;
    for (size_t i = 0; i < traces.size(); i++) {
        for (size_t j = 0; j < traces[i].size(); j++) {
            traces[i][j].z = heights.interpolate(traces[i][j]) + .5f;
        }
    }
    return traces;
}

void Grid::thermalErosion()
{

}

void Grid::windErosion()
{

}

void Grid::randomFaultTerrainGeneration(int numberOfFaults, int maxNumberOfSubpointsInFaults, float faultHeight)
{
    std::vector<BSpline> faults;
    for (int i = 0; i < numberOfFaults; i++) {
        Vector3 firstPoint, lastPoint;
        int numberOfSubpoints = random_gen::generate(0, maxNumberOfSubpointsInFaults + 1);

        // Force first and last points to be on the borders of the grid (maybe not the mendatory?)
        if (random_gen::generate() < .5f) firstPoint = Vector3(random_gen::generate(0, getSizeX()), std::round(random_gen::generate()) * getSizeY());
        else                            firstPoint = Vector3(std::round(random_gen::generate()) * getSizeX(), random_gen::generate(0, getSizeY()));
        if (random_gen::generate() < .5f) lastPoint = Vector3(random_gen::generate(0, getSizeX()), std::round(random_gen::generate()) * getSizeY());
        else                            lastPoint = Vector3(std::round(random_gen::generate()) * getSizeX(), random_gen::generate(0, getSizeY()));

        std::vector<Vector3> points;
        points.push_back(firstPoint);
        for (int iSub = 0; iSub < numberOfSubpoints; iSub++)
            points.push_back(Vector3(random_gen::generate(0, getSizeX()), random_gen::generate(0, getSizeY())));
        points.push_back(lastPoint);
        faults.push_back(BSpline(points));
    }

    Matrix3<float> faultsImpact(getSizeX(), getSizeY());
    auto start = std::chrono::system_clock::now();
    for (int x = 0; x < getSizeX(); x++) {
        for (int y = 0; y < getSizeY(); y++) {
            Vector3 pos(x, y);
            float totalInfluence = 0;
            for (size_t i = 0; i < faults.size(); i++) {
                float coef = 1.f; // / (float)(i + 1);
                float distanceToBorder = std::min(getSizeX() - x, std::min(x, std::min(getSizeY() - y, y)));
                totalInfluence += std::max(0.f, interpolation::fault_distance(faults[i].estimateDistanceFrom(pos), 2.f) - interpolation::fault_distance(distanceToBorder, 20.f)) * coef;
            }
            faultsImpact.at(pos) = totalInfluence / (float)numberOfFaults;
        }
    }
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count() << "ms" << std::endl;
    this->heights += (faultsImpact.normalized() * faultHeight);
}
void Grid::fromVoxelGrid(VoxelGrid &voxelGrid) {
    for (int x = 0; x < voxelGrid.getSizeX(); x++) {
        for (int y = 0; y < voxelGrid.getSizeY(); y++) {
            this->heights.at(x, y) = voxelGrid.getHeight(x, y);
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
    toFloatData = (this->heights/newHeight).data;
    for (size_t i = 0; i < this->heights.size(); i++) {
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
