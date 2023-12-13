#include "DataStructure/Voxel.h"

#include <sstream>
#include "DataStructure/Matrix3.h"

VoxelDataFile::VoxelDataFile()
{}

VoxelDataFile::VoxelDataFile(int w, int h, int d, const GridF&dVec)
    : width(w), height(h), depth(d), data(dVec) {}

void VoxelDataFile::write(const std::string &filename) {
    std::ofstream outFile(filename, std::ios::binary);
    if (outFile.is_open()) {
        outFile.write(reinterpret_cast<const char*>(&width), sizeof(width));
        outFile.write(reinterpret_cast<const char*>(&height), sizeof(height));
        outFile.write(reinterpret_cast<const char*>(&depth), sizeof(depth));

        size_t dataSize = data.size();
        outFile.write(reinterpret_cast<const char*>(&dataSize), sizeof(dataSize));
        outFile.write(reinterpret_cast<const char*>(data.data.data()), dataSize * sizeof(float));

        outFile.close();
        std::cout << "Data written to file successfully." << std::endl;
    } else {
        std::cerr << "Unable to open file for writing." << std::endl;
    }
}

void VoxelDataFile::load(const std::string &filename) {
    std::ifstream inFile(filename);
    if (inFile.is_open()) {
        std::string line;
        if (std::getline(inFile, line)) {
            std::istringstream iss(line);
            int w, h, d, dump;
            iss >> w >> h >> d >> dump;

            // Check for the old text format
            if (iss.eof() && w > 0 && h > 0 && d > 0) {
                // Rewind and load using the old text format
                inFile.seekg(0, std::ios::beg);
                loadFromFileOld(filename);
                return;
            }
        }

        // Rewind and check for the new binary format
        inFile.seekg(0, std::ios::beg);
        int potentialWidth, potentialHeight, potentialDepth;
        inFile.read(reinterpret_cast<char*>(&potentialWidth), sizeof(potentialWidth));
        inFile.read(reinterpret_cast<char*>(&potentialDepth), sizeof(potentialDepth));
        inFile.read(reinterpret_cast<char*>(&potentialHeight), sizeof(potentialHeight));

        if (potentialWidth > 0 && potentialHeight > 0 && potentialDepth > 0) {
            // Assume it's the new binary format
            inFile.seekg(0, std::ios::beg);
            loadFromFileBinary(filename);
        } else {
            std::cerr << "Unable to determine the file format." << std::endl;
        }

        inFile.close();
    } else {
        std::cerr << "Unable to open file for reading." << std::endl;
    }
}

void VoxelDataFile::loadFromFileBinary(const std::string &filename) {
    std::ifstream inFile(filename, std::ios::binary);
    if (inFile.is_open()) {
        inFile.read(reinterpret_cast<char*>(&width), sizeof(width));
        inFile.read(reinterpret_cast<char*>(&depth), sizeof(depth));
        inFile.read(reinterpret_cast<char*>(&height), sizeof(height));

        size_t dataSize;
        inFile.read(reinterpret_cast<char*>(&dataSize), sizeof(dataSize));
        data.data.resize(dataSize);
        inFile.read(reinterpret_cast<char*>(data.data.data()), dataSize * sizeof(float));
        data.sizeX = width;
        data.sizeY = depth;
        data.sizeZ = height;

        inFile.close();
        std::cout << "Data loaded from file successfully." << std::endl;
    } else {
        std::cerr << "Unable to open file for reading." << std::endl;
    }
}

void VoxelDataFile::loadFromFileOld(const std::string &filename) {
    std::ifstream in;
    in.open(filename);
    if (in.fail()) {
        std::cerr << "Unable to open file " << filename << "..." << std::endl;
        return;
    }
    int _x, _y, _z;
    int _chunkSize;

    std::string firstLine;
    std::getline(in, firstLine);  // Read the entire first line

    std::stringstream ss(firstLine);
    ss >> _x >> _y >> _z;


    Vector3 finalSize = Vector3(_x, _y, _z);
    GridF _cachedVoxelValues = GridF(_x, _y, _z, 0.f);

    if (ss >> _chunkSize) { // Compatibility with previous terrain saving system
        int chunkSize = _chunkSize;
        float map_val;
        int iChunk = 0;
        for (int xChunk = 0; xChunk < std::ceil(_x / (float)chunkSize); xChunk++) {
            for (int yChunk = 0; yChunk < std::ceil(_y / (float)chunkSize); yChunk++) {
                Vector3 offset(xChunk * chunkSize, yChunk * chunkSize, 0.f);
                for (int x = 0; x < chunkSize; x++) {
                    for (int y = 0; y < chunkSize; y++) {
                        for (int z = 0; z < _z; z++) {
                            in >> map_val;
                            _cachedVoxelValues.at(Vector3(x, y, z) + offset) = map_val;
                        }
                    }
                }
            }
        }
    } else {
        std::cout << "Retrieving..." << std::endl;
        for (size_t i = 0; i < _cachedVoxelValues.size(); i++)
            in >> _cachedVoxelValues[i];
        std::cout << _cachedVoxelValues << std::endl;
    }
    _cachedVoxelValues.iterateParallel([&](size_t i) {
        if (_cachedVoxelValues[i] < -1.f)
            _cachedVoxelValues[i] = -1.f;
    });
    this->width = _x;
    this->depth = _y;
    this->height = _z;
    this->data = _cachedVoxelValues;
}
