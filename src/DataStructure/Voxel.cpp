#include "DataStructure/Voxel.h"

#include <sstream>
#include "DataStructure/Matrix3.h"

VoxelDataFile::VoxelDataFile()
{}

VoxelDataFile::VoxelDataFile(const GridF&dVec)
    : data(dVec) {}

void VoxelDataFile::write(const std::string &filename) {
    std::ofstream outFile(filename, std::ios::binary);
    if (outFile.is_open()) {
        int width = data.sizeX, depth = data.sizeY, height = data.sizeZ;
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
            try {
                inFile.seekg(0, std::ios::beg);
                loadFromFileBinary(filename);
            } catch(std::exception e) {
                inFile.seekg(0, std::ios::beg);
                loadFromFileOld(filename);
            }
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
        int width, depth, height;
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

    GridF _cachedVoxelValues = GridF(_x, _y, _z, 0.f);

    if (ss >> _chunkSize) { // Compatibility with previous terrain saving system
        int chunkSize = _chunkSize;
        float map_val;
//        int iChunk = 0;
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
    }
    _cachedVoxelValues.iterateParallel([&](size_t i) {
        if (_cachedVoxelValues[i] < -1.f)
            _cachedVoxelValues[i] = -1.f;
    });
    this->data = _cachedVoxelValues;
}








VectorFieldDataFile::VectorFieldDataFile()
{}

VectorFieldDataFile::VectorFieldDataFile(const GridV3& dVec)
    : data(dVec) {}

void VectorFieldDataFile::write(const std::string &filename) {
    std::ofstream outFile(filename, std::ios::binary);
    if (outFile.is_open()) {
        int width = data.sizeX, depth = data.sizeY, height = data.sizeZ;
        outFile.write(reinterpret_cast<const char*>(&width), sizeof(width));
        outFile.write(reinterpret_cast<const char*>(&depth), sizeof(depth));
        outFile.write(reinterpret_cast<const char*>(&height), sizeof(height));

        size_t dataSize = data.size();
        outFile.write(reinterpret_cast<const char*>(&dataSize), sizeof(dataSize));
        // Write each Vector3 as raw binary data
        for (const auto& vec : data.data) {
            outFile.write(reinterpret_cast<const char*>(&vec), sizeof(vec));
        }

        outFile.close();
//        std::cout << "Data written to file successfully." << std::endl;
    } else {
        std::cerr << "Unable to open file " << filename << " for writing." << std::endl;
    }
}

void VectorFieldDataFile::load(const std::string &filename) {
    loadFromFileBinary(filename);
}

void VectorFieldDataFile::loadFromFileBinary(const std::string &filename) {
    std::ifstream inFile(filename, std::ios::binary);
    if (inFile.is_open()) {
        int width, depth, height;
        inFile.read(reinterpret_cast<char*>(&width), sizeof(width));
        inFile.read(reinterpret_cast<char*>(&depth), sizeof(depth));
        inFile.read(reinterpret_cast<char*>(&height), sizeof(height));

        size_t dataSize;
        inFile.read(reinterpret_cast<char*>(&dataSize), sizeof(dataSize));
        data.data.resize(dataSize);
        // Read each Vector3 from raw binary data
        for (size_t i = 0; i < dataSize; ++i) {
            auto& vec = data.data[i];
            inFile.read(reinterpret_cast<char*>(&vec), sizeof(vec));/*
            inFile.read(reinterpret_cast<char*>(&(vec.x)), sizeof(vec.x));
            inFile.read(reinterpret_cast<char*>(&(vec.y)), sizeof(vec.y));
            inFile.read(reinterpret_cast<char*>(&(vec.z)), sizeof(vec.z));
            std::cout << vec << " " << std::flush;*/
//            inFile.read(reinterpret_cast<const char*>(&vec.valid), sizeof(vec.valid));
        }
        data.sizeX = width;
        data.sizeY = depth;
        data.sizeZ = height;

        inFile.close();
        std::cout << "Data loaded from file successfully." << std::endl;
    } else {
        std::cerr << "Unable to open file for reading." << std::endl;
    }
}
