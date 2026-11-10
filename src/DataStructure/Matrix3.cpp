#include "DataStructure/Matrix3.h"
//#include "Utils/stb_image.h"
#include "DataStructure/Image.h"
#include "Utils/Skeletonize.h"
#include <queue>
#include <sstream>
#include <cstring>
#include <locale>

template <>
Matrix3<Vector3> Matrix3<Vector3>::curl(float radius) const {
    if (!std::isfinite(radius) || radius <= 0)
        throw std::invalid_argument("Curl radius must be positive and finite");
    Matrix3<Vector3> returningGrid(this->sizeX, this->sizeY, this->sizeZ);
//    float radius = 1;
    iterateParallel([&] (size_t x, size_t y, size_t z) {
        // Central finite difference
        Vector3 dX = at(int(x) + radius, y, z) - at(int(x) - radius, y, z);
        Vector3 dY = at(x, int(y) + radius, z) - at(x, int(y) - radius, z);
        Vector3 dZ = at(x, y, int(z) + radius) - at(x, y, int(z) - radius);

        Vector3 curl(dY.z() - dZ.y(), dZ.x() - dX.z(), dX.y() - dY.x());
        returningGrid(x, y, z) = curl / (radius * 2.f);
    });
    return returningGrid;
}
template <>
inline Matrix3<Vector3> Matrix3<Vector3>::rot() const {
    return this->curl();
}

template <>
Matrix3<int> Matrix3<int>::skeletonize() const
{
    Matrix3<int> self = ((Matrix3<float>)*this).binarize(0.5f);
    Matrix3<int> initial = *this;
    skeleton_tracer_t* skel = new skeleton_tracer_t();
    skel->W = self.sizeX; // width of image
    skel->H = self.sizeY; // height of image

    // allocate the input image
    unsigned char* data = (unsigned char*)malloc(sizeof(unsigned char)*skel->W*skel->H); //new uchar(self.sizeX * self.sizeY);
    self.iterateParallel([&] (size_t i) {
        data[i] = (unsigned char)(self[i]);
    });
    /*for (size_t i = 0; i < self.size(); i++)
        data[i] = (unsigned char)(self[i]);*/
    skel->im = data;

    skel->thinning_zs(); // perform raster thinning

    // run the algorithm
    skeleton_tracer_t::polyline_t* p = (skeleton_tracer_t::polyline_t*)skel->trace_skeleton(0, 0, skel->W, skel->H, 0);

    self.reset(); // Back to all 0's
    // print out points in every polyline
    skeleton_tracer_t::polyline_t* it = p; //iterator
    while(it){
      skeleton_tracer_t::point_t* jt = it->head;
      Vector3 prevPos(false);
      while(jt){
          Vector3 newPos(jt->x, jt->y);
          // Clearly not the best way to do this!
          if (prevPos.isValid()/* && initial.at(newPos) != 0*/) {
              Vector3 move = newPos - prevPos;
              float length = move.length();
              for (int i = 0; i < length; i++)
                  self.at(prevPos + move * clamp(i / length, 0.f, 1.f)) = 1;
          }
          prevPos = newPos;
        jt = jt->next;
      }
      it = it->next;
    }
    free(skel->im);
    skel->destroy_polylines(p);
    skel->destroy_rects();
    delete skel;
    return self;
}

template <>
std::vector<Polyline> Matrix3<int>::skeletonizeToBSplines() const
{
    Matrix3<int> initial = ((Matrix3<float>)*this).binarize(0.5f);
    skeleton_tracer_t* skel = new skeleton_tracer_t();
    skel->W = this->sizeX; // width of image
    skel->H = this->sizeY; // height of image

    // allocate the input image
    unsigned char* data = (unsigned char*)malloc(sizeof(unsigned char)*skel->W*skel->H); //new uchar(self.sizeX * self.sizeY);
    iterateParallel([&] (size_t i) {
        data[i] = (unsigned char)(initial(i));
    });
    /*for (size_t i = 0; i < self.size(); i++)
        data[i] = (unsigned char)(self[i]);*/
    skel->im = data;

    skel->thinning_zs(); // perform raster thinning

    // run the algorithm
    skeleton_tracer_t::polyline_t* p = (skeleton_tracer_t::polyline_t*)skel->trace_skeleton(0, 0, skel->W, skel->H, 0);


    std::vector<std::vector<Vector3>> splines;
    skeleton_tracer_t::polyline_t* it = p; //iterator
    while(it){
      skeleton_tracer_t::point_t* jt = it->head;
        std::vector<Vector3> spline;
      while(jt){
          ;
          spline.push_back(Vector3(jt->x, jt->y));
          jt = jt->next;
      }
      it = it->next;
      splines.push_back(spline);
    }
    free(skel->im);
    skel->destroy_polylines(p);
    skel->destroy_rects();
    delete skel;
//    return splines;

    // Try at best to merge the curves
    float limitDistance = 20.f;
    float sqrLim = limitDistance * limitDistance;
    bool atLeastOneMerging = true;
    while (atLeastOneMerging) {
        atLeastOneMerging = false;
        std::vector<std::vector<Vector3>> merged;
        for (int i = int(splines.size()) - 1; i >= 0; i--) {
            auto& spline = splines[i];
            bool isMerged = false;
            for (auto& merge : merged) {
                // Try back to front, front to back, back to back and front to front
                auto frontSpline = spline.front();
                auto backSpline = spline.back();
                auto frontMerge = merge.front();
                auto backMerge = merge.back();
                // back to front:
                if ((frontMerge - backSpline).norm2() < sqrLim) {
                    merge.insert(merge.begin(), spline.begin(), spline.end());
                    isMerged = true;
                    break;
                }
                // front to back:
                else if ((backMerge - frontSpline).norm2() < sqrLim) {
                    merge.insert(merge.end(), spline.begin(), spline.end());
                    isMerged = true;
                    break;
                }
                // back to back:
                else if ((backMerge - backSpline).norm2() < sqrLim) {
                    std::reverse(spline.begin(), spline.end());
                    merge.insert(merge.end(), spline.begin(), spline.end());
                    isMerged = true;
                    break;
                }
                // front to front:
                else if ((frontMerge - frontSpline).norm2() < sqrLim) {
                    std::reverse(spline.begin(), spline.end());
                    merge.insert(merge.begin(), spline.begin(), spline.end());
                    isMerged = true;
                    break;
                }
            }
            if (!isMerged)
                merged.push_back(spline);
            else
                atLeastOneMerging = true;
        }
        splines = merged;
    }
    std::vector<Polyline> allSplines;
    allSplines.reserve(splines.size());
    for (const auto& spline : splines) { allSplines.push_back(spline); }
    return allSplines;
}
template <>
Matrix3<Vector3> Matrix3<float>::fillWithBSplines(std::vector<Polyline> splines) const {
    Matrix3<Vector3> newImage(this->getDimensions());
    newImage.raiseErrorOnBadCoord = false;
    int pointsize = 5;
    int i = 0;
    for (auto& spline : splines) {
        for (auto& point : spline) {
            for (int x = point.x() - pointsize; x < point.x() + pointsize; x++) {
                for (int y = point.y() - pointsize; y < point.y() + pointsize; y++) {
                    newImage.at(x, y) = HSVtoRGB((float)i / splines.size(), 1.0f, 0.5f);
                }
            }
        }
        for (auto& point : spline.resamplePoints(10000)) {
            // newImage.at(point) = HSVtoRGB((float)i / 10000.f, 1.0f, 1.0f);
            newImage.at(point) = HSVtoRGB((float)i / splines.size(), 1.0f, 1.0f);
        }
        i++;
    }
    newImage.raiseErrorOnBadCoord = true;
    return newImage;
}

template <>
Matrix3<int> Matrix3<int>::computeConnectedComponents(bool use4Connect) const
{
    Matrix3<int> result(getDimensions());
    int label = 0;
    std::vector<Vector3i> stack;
    iterate([&](int x, int y, int z) {
        Vector3i p(x, y, z);
        if (!unsafe(p) || result.unsafe(p))
            return;
        ++label;
        result.unsafe(p) = label;
        stack.push_back(p);
        while (!stack.empty()) {
            auto current = stack.back();
            stack.pop_back();
            for (int dz = -1; dz <= 1; ++dz)
                for (int dy = -1; dy <= 1; ++dy)
                    for (int dx = -1; dx <= 1; ++dx) {
                        int distance = std::abs(dx) + std::abs(dy) + std::abs(dz);
                        if (distance == 0 || (use4Connect && distance != 1))
                            continue;
                        Vector3i q(current.x() + dx, current.y() + dy, current.z() + dz);
                        if (checkCoord(q) && unsafe(q) && !result.unsafe(q)) {
                            result.unsafe(q) = label;
                            stack.push_back(q);
                        }
                    }
        }
    });
    return result;
}

template <>
Matrix3<int> Matrix3<int>::findContour(bool use2D) const
{
    auto eroded = this->erode(use2D);
    return *this - eroded;
}

template <>
std::vector<Contour> Matrix3<int>::findContoursAsCurves() const
{
    std::vector<Contour> curves;
    auto grid = this->resize(sizeX * 2.f, sizeY * 2.f, RESIZE_MODE::MAX_VAL).findContour(true);
//    std::cout << "Grid : " << grid << "\n" << this->resize(sizeX * .5f, sizeY * .5f, RESIZE_MODE::MAX_VAL).findContour(true).displayAsPlot() << std::endl;
    grid.raiseErrorOnBadCoord = false;

    while (grid.max() != 0) {
        std::vector<Vector3> contour;
        // Find the starting point (the first 1 encountered in the grid)
        Vector3i start(false);

        grid.iterate([&] (size_t i) {
            if (!start.isValid() && grid[i] == 1) {
                start = grid.getCoordAsVector3(i);
                return;
            }
        });

        // If no point was found, return the empty contour
        if (!start.isValid()) return curves;

        // Start from the initial point and do a simple DFS to find the contour
        std::vector<Vector3i> q;
        q.push_back(start);
        grid(start) = 0;  // Mark as visited

        while (!q.empty()) {
            Vector3i current = q.back();
            q.pop_back();
            // Add to contour
            contour.push_back(current);

            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    if (!(dx == 0 || dy == 0)) continue;
                    Vector3i nextPos(current.x() + dx, current.y() + dy);
                    if (grid(nextPos)) {
                        q.push_back(nextPos);
                        grid(nextPos) = 0;
                    }
                }
            }
        }

        for (auto &p : contour)
            p *= .5f;
        curves.push_back(
            Polyline(contour) /*.resamplePoints().simplifyByRamerDouglasPeucker(5.f).resamplePoints()*/);
    }
    return curves;
}

template <> Matrix3<float> Matrix3<Vector3>::divergence() const {
    auto self = *this;
    self.raiseErrorOnBadCoord = false;
    Matrix3<float> returningGrid(this->sizeX, this->sizeY, this->sizeZ);
    iterateParallel([&](int x, int y, int z) {
        returningGrid.at(x, y, z) = ((self.at(x + 1, y, z) - self.at(x - 1, y, z)).x() +
                                     (self.at(x, y + 1, z) - self.at(x, y - 1, z)).y() +
                                     (self.at(x, y, z + 1) - self.at(x, y, z - 1)).z()) *
                                    .5f;
    });
    return returningGrid;
}

template <> Vector3 Matrix3<Vector3>::gradient(const Vector3 &position) const {
    auto self = *this;
    self.raiseErrorOnBadCoord = false;
    Vector3 flooredPos = position.floor();
    Vector3 offset = position - flooredPos;
    return Vector3(
        self.at(flooredPos + Vector3(1, 0, 0)).x() * (1 - offset.x()) + self.at(flooredPos).x() * offset.x(),
        self.at(flooredPos + Vector3(0, 1, 0)).y() * (1 - offset.y()) + self.at(flooredPos).y() * offset.y(),
        self.at(flooredPos + Vector3(0, 0, 1)).z() * (1 - offset.z()) + self.at(flooredPos).z() * offset.z());
}

template <> inline Vector3 Matrix3<Vector3>::gradient(float posX, float posY, float posZ) const {
    return gradient(Vector3(posX, posY, posZ));
}

template <> Matrix3<Vector3> Matrix3<Vector3>::gradient() const {
    auto self = *this;
    self.raiseErrorOnBadCoord = false;
    Matrix3<Vector3> returningGrid(this->sizeX, this->sizeY, this->sizeZ);
    iterateParallel([&](int x, int y, int z) {
        returningGrid.at(x, y, z) = Vector3((self.at(x + 1, y, z) - self.at(x - 1, y, z)).x() * .5f,
                                            (self.at(x, y + 1, z) - self.at(x, y - 1, z)).y() * .5f,
                                            (self.at(x, y, z + 1) - self.at(x, y, z - 1)).z() * .5f);
    });
    return returningGrid;
}

template <> Matrix3<Vector3> Matrix3<Vector3>::random(size_t sizeX, size_t sizeY, size_t sizeZ) {
    Matrix3<Vector3> result(sizeX, sizeY, sizeZ);
    for (auto &value : result.data)
        value = Vector3::random();
    return result;
}

template <> Matrix3<Vector3> Matrix3<Vector3>::fromImageRGB(const std::string &filename) {
    auto image = Image::readFromFile(filename);
    if (image.isColor())
        return std::move(image.colorImage);
    return image.getColorImage();
}

template <> Matrix3<float> Matrix3<float>::fromImageBW(const std::string &filename) {
    int width = 0, height = 0, channels = 0;
    if (stbi_is_hdr(filename.c_str()))
        return Image::readFromFile(filename).getBwImage();
    if (stbi_is_16_bit(filename.c_str())) {
        std::unique_ptr<stbi_us, decltype(&stbi_image_free)> pixels(
            stbi_load_16(filename.c_str(), &width, &height, &channels, 1), stbi_image_free);
        if (!pixels)
            throw std::runtime_error("Cannot decode grayscale image: " + filename);
        GridF result(width, height);
        result.iterateParallel([&](size_t i) { result.data[i] = pixels.get()[i] / 65535.f; });
        return result;
    }
    std::unique_ptr<stbi_uc, decltype(&stbi_image_free)> pixels(
        stbi_load(filename.c_str(), &width, &height, &channels, 1), stbi_image_free);
    if (!pixels)
        throw std::runtime_error("Cannot decode grayscale image: " + filename);
    GridF result(width, height);
    result.iterateParallel([&](size_t i) { result.data[i] = pixels.get()[i] / 255.f; });
    return result;
}
// template <class T>
Matrix3<float> operator-(const float a, const Matrix3<float> &b) {
    Matrix3<float> res = b;
    res.iterateParallel([&](size_t i) { res[i] = a - res[i]; });
    /*for (size_t i = 0; i < res.size(); i++)
        res[i] = a - res[i];*/
    return res;
}
// template <class T>
Matrix3<float> operator+(const float a, const Matrix3<float> &b) {
    return b + a;
    /*Matrix3<float> res = b;
    for (size_t i = 0; i < res.size(); i++)
        res[i] = a + res[i];
    return res;*/
}

namespace {
constexpr size_t gridHeaderSize = 3 * sizeof(int) + sizeof(size_t);
void readHeader(const std::string &bytes, int &x, int &y, int &z, size_t &count) {
    if (bytes.size() < gridHeaderSize)
        throw std::invalid_argument("Truncated grid header");
    std::memcpy(&x, bytes.data(), sizeof(int));
    std::memcpy(&y, bytes.data() + sizeof(int), sizeof(int));
    std::memcpy(&z, bytes.data() + 2 * sizeof(int), sizeof(int));
    std::memcpy(&count, bytes.data() + 3 * sizeof(int), sizeof(size_t));
    if (x < 0 || y < 0 || z < 0)
        throw std::invalid_argument("Negative grid dimensions");
    if (count != grid_detail::count<float>(x, y, z))
        throw std::invalid_argument("Grid payload count mismatch");
}
std::string gridHeader(const Vector3i &dims, size_t count) {
    std::string bytes(gridHeaderSize, '\0');
    int x = dims.x(), y = dims.y(), z = dims.z();
    std::memcpy(bytes.data(), &x, sizeof(x));
    std::memcpy(bytes.data() + sizeof(x), &y, sizeof(y));
    std::memcpy(bytes.data() + 2 * sizeof(x), &z, sizeof(z));
    std::memcpy(bytes.data() + 3 * sizeof(x), &count, sizeof(count));
    return bytes;
}
void payloadLength(const std::string &bytes, size_t count, size_t itemSize) {
    if (count > (std::numeric_limits<size_t>::max() - gridHeaderSize) / itemSize ||
        bytes.size() != gridHeaderSize + count * itemSize)
        throw std::invalid_argument("Truncated or oversized grid payload");
}
void textEnd(std::istringstream &in) {
    in >> std::ws;
    if (!in.eof())
        throw std::invalid_argument("Unexpected trailing grid data");
}
Vector3i textHeader(std::istringstream &in) {
    int x, y, z;
    if (!(in >> x >> y >> z) || x < 0 || y < 0 || z < 0)
        throw std::invalid_argument("Invalid text grid header");
    grid_detail::count<float>(x, y, z);
    return Vector3i(x, y, z);
}
} // namespace
std::string stringifyGridF(const GridF &data, bool binaryMode) {
    if (data.size() != grid_detail::count<float>(data.sizeX, data.sizeY, data.sizeZ))
        throw std::invalid_argument("Inconsistent grid storage");
    if (binaryMode) {
        auto bytes = gridHeader(data.getDimensions(), data.size());
        if (!data.empty())
            bytes.append(reinterpret_cast<const char *>(data.data.data()), data.size() * sizeof(float));
        return bytes;
    }
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::setprecision(std::numeric_limits<float>::max_digits10);
    out << data.sizeX << ' ' << data.sizeY << ' ' << data.sizeZ;
    for (float v : data)
        out << ' ' << v;
    return out.str();
}
std::string stringifyGridV3(const GridV3 &data, bool binaryMode) {
    if (data.size() != grid_detail::count<Vector3>(data.sizeX, data.sizeY, data.sizeZ))
        throw std::invalid_argument("Inconsistent grid storage");
    if (binaryMode) {
        // Preserve existing native-format files, but never emit object padding.
        auto bytes = gridHeader(data.getDimensions(), data.size());
        bytes.resize(gridHeaderSize + data.size() * sizeof(Vector3), 0);
        for (size_t i = 0; i < data.size(); ++i) {
            auto dst = bytes.data() + gridHeaderSize + i * sizeof(Vector3);
            std::memcpy(dst + offsetof(Vector3, v), data[i].data(), 3 * sizeof(float));
            bool valid = data[i].valid;
            std::memcpy(dst + offsetof(Vector3, valid), &valid, sizeof(bool));
        }
        return bytes;
    }
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::setprecision(std::numeric_limits<float>::max_digits10);
    out << data.sizeX << ' ' << data.sizeY << ' ' << data.sizeZ;
    for (auto v : data)
        out << ' ' << v.x() << ' ' << v.y() << ' ' << v.z();
    return out.str();
}
GridF loadGridF(const std::string &bytes, bool binaryMode) {
    if (binaryMode) {
        int x, y, z;
        size_t count;
        readHeader(bytes, x, y, z, count);
        payloadLength(bytes, count, sizeof(float));
        GridF result(x, y, z);
        if (count)
            std::memcpy(result.data.data(), bytes.data() + gridHeaderSize, count * sizeof(float));
        return result;
    }
    std::istringstream in(bytes);
    in.imbue(std::locale::classic());
    auto dims = textHeader(in);
    // At least one textual character is needed for every value; validate before allocation.
    if (grid_detail::count<float>(dims.x(), dims.y(), dims.z()) > bytes.size())
        throw std::invalid_argument("Truncated text grid");
    GridF result(dims);
    for (auto &v : result)
        if (!(in >> v))
            throw std::invalid_argument("Invalid text sample");
    textEnd(in);
    return result;
}
GridV3 loadGridV3(const std::string &bytes, bool binaryMode) {
    if (binaryMode) {
        int x, y, z;
        size_t count;
        readHeader(bytes, x, y, z, count);
        payloadLength(bytes, count, sizeof(Vector3));
        GridV3 result(x, y, z);
        for (size_t i = 0; i < count; ++i) {
            auto src = bytes.data() + gridHeaderSize + i * sizeof(Vector3);
            // Deserialize arithmetic components, never copy untrusted bytes into a bool.
            const unsigned char flag = static_cast<unsigned char>(src[offsetof(Vector3, valid)]);
            if (flag > 1)
                throw std::invalid_argument("Invalid vector validity flag");
            std::memcpy(result[i].data(), src + offsetof(Vector3, v), 3 * sizeof(float));
            result[i].valid = flag != 0;
        }
        return result;
    }
    std::istringstream in(bytes);
    in.imbue(std::locale::classic());
    auto dims = textHeader(in);
    if (grid_detail::count<Vector3>(dims.x(), dims.y(), dims.z()) > bytes.size() / 3)
        throw std::invalid_argument("Truncated text vector grid");
    GridV3 result(dims);
    for (auto &v : result)
        if (!(in >> v.x() >> v.y() >> v.z()))
            throw std::invalid_argument("Invalid vector sample");
    textEnd(in);
    return result;
}
