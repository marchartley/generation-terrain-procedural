#include "DataStructure/Voxel.h"
#include <sstream>
#include <limits>
#include <filesystem>
#include <cstdint>
#include <QSaveFile>

namespace {
std::string readFile(const std::string &name) {
    std::ifstream in(name, std::ios::binary);
    if (!in)
        throw std::runtime_error("Cannot open grid file: " + name);
    std::ostringstream bytes;
    bytes << in.rdbuf();
    if (in.bad())
        throw std::runtime_error("Grid file read failed");
    return bytes.str();
}
void writeFile(const std::string &name, const std::string &bytes) {
    QSaveFile out(QString::fromStdString(name));
    if (!out.open(QIODevice::WriteOnly))
        throw std::runtime_error("Cannot open grid output: " + name);
    if (out.write(bytes.data(), qint64(bytes.size())) != qint64(bytes.size()) || !out.commit())
        throw std::runtime_error("Atomic grid output failed: " + name);
}
} // namespace
VoxelDataFile::VoxelDataFile() = default;
VoxelDataFile::VoxelDataFile(const GridF &values) : data(values) {}
void VoxelDataFile::write(const std::string &name) { writeFile(name, stringifyGridF(data)); }
void VoxelDataFile::load(const std::string &name) {
    auto bytes = readFile(name);
    std::istringstream header(bytes.substr(0, bytes.find('\n')));
    int x, y, z;
    if ((header >> x >> y >> z) && x >= 0 && y >= 0 && z >= 0)
        loadFromFileOld(name);
    else
        data = loadGridF(bytes);
}
void VoxelDataFile::loadFromFileBinary(const std::string &name) {
    auto replacement = loadGridF(readFile(name));
    data = std::move(replacement);
}
void VoxelDataFile::loadFromFileOld(const std::string &name) {
    std::ifstream in(name);
    std::string line;
    std::getline(in, line);
    std::istringstream header(line);
    int x, y, z, chunk = 0;
    if (!(header >> x >> y >> z) || x < 0 || y < 0 || z < 0)
        throw std::invalid_argument("Invalid legacy voxel header");
    const auto count = grid_detail::count<float>(x, y, z);
    const auto fileBytes = std::filesystem::file_size(name);
    if (count > fileBytes)
        throw std::invalid_argument("Legacy voxel payload is too short for its dimensions");
    header >> std::ws;
    const bool chunked = !header.eof();
    if (chunked && !(header >> chunk))
        throw std::invalid_argument("Invalid legacy chunk header");
    header >> std::ws;
    if (!header.eof())
        throw std::invalid_argument("Trailing legacy header fields");
    if (chunked && chunk <= 0)
        throw std::invalid_argument("Invalid legacy chunk size");
    if (chunked && count) {
        uintmax_t required = 1;
        for (uintmax_t factor : {(uintmax_t(x) + chunk - 1) / chunk, (uintmax_t(y) + chunk - 1) / chunk,
                                 uintmax_t(chunk), uintmax_t(chunk), uintmax_t(z)}) {
            if (factor && required > fileBytes / factor)
                throw std::invalid_argument("Truncated padded legacy chunks");
            required *= factor;
        }
    }
    GridF replacement(x, y, z);
    if (chunked) {
        if (chunk <= 0)
            throw std::invalid_argument("Invalid legacy chunk size");
        for (int64_t cx = 0; cx < x; cx += chunk)
            for (int64_t cy = 0; cy < y; cy += chunk)
                for (int64_t dx = 0; dx < chunk; ++dx)
                    for (int64_t dy = 0; dy < chunk; ++dy)
                        for (int zz = 0; zz < z; ++zz) {
                            float value;
                            if (!(in >> value))
                                throw std::invalid_argument("Truncated legacy chunks");
                            if (cx + dx < x && cy + dy < y)
                                replacement.unsafe(cx + dx, cy + dy, zz) = std::max(-1.f, value);
                        }
    } else {
        for (auto &value : replacement) {
            if (!(in >> value))
                throw std::invalid_argument("Truncated legacy voxel data");
            value = std::max(-1.f, value);
        }
    }
    in >> std::ws;
    if (!in.eof())
        throw std::invalid_argument("Trailing legacy voxel data");
    data = std::move(replacement);
}
VectorFieldDataFile::VectorFieldDataFile() = default;
VectorFieldDataFile::VectorFieldDataFile(const GridV3 &values) : data(values) {}
void VectorFieldDataFile::write(const std::string &name) { writeFile(name, stringifyGridV3(data)); }
void VectorFieldDataFile::load(const std::string &name) { loadFromFileBinary(name); }
void VectorFieldDataFile::loadFromFileBinary(const std::string &name) {
    auto replacement = loadGridV3(readFile(name));
    data = std::move(replacement);
}
