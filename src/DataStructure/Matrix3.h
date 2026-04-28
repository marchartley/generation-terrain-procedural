#ifndef MATRIX3_H
#define MATRIX3_H

#include <iostream>
#include <vector>
#include <tuple>
#include <memory>
#include <string>
#include <math.h>
#include <iomanip>
#include <complex>
#include "DataStructure/Vector3.h"
#include "Curves/CatmullRomSpline.h"
#include "Curves/ShapeCurve.h"
#include "Utils/Collisions.h"
#include "Utils/Utils.h"

#include "Utils/FastNoiseLit.h"

enum RESIZE_MODE {
    LINEAR = 0,
    MAX_VAL = 1,
    MIN_VAL = 2,
    FILL_WITH_DEFAULT = 3,
    NEAREST = 4
};

enum CONVOLUTION_BORDERS {
    ZERO_PAD = 0,
    REPEAT = 1,
    MIRROR = 2,
    WRAPPING = 3,
    IGNORED = 4,
    COPY = 5
};

enum RETURN_VALUE_ON_OUTSIDE {
    DEFAULT_VALUE = 0,
    MIRROR_VALUE = 1,
    WRAPPED_VALUE = 2,
    REPEAT_VALUE = 3
};

enum NORMALIZE_METHOD {
    NORMALIZE_MINMAX,
    NORMALIZE_Z_SCORE,
    NORMALIZE_SOFTMAX
};

// Warning : don't use bool type...
// This class is based on a std::vector, which has specifications on bools
// the [] operator won't work ... Use int or short int instead
template <class T>
class Matrix3
{
public:
    Matrix3();
    Matrix3(size_t sizeX, size_t sizeY, size_t sizeZ = 1, T initValue = T());
    // Matrix3(const Vector3& size, T initValue = T());
    Matrix3(const Vector3i& size, T initValue = T());
    Matrix3(const std::vector<std::vector<std::vector<T>>>& data);
    Matrix3(const std::vector<std::vector<T>>& data);
    Matrix3(const std::vector<T>& data, size_t sizeX, size_t sizeY, int sizeZ = -1);

    const T& at(int i, int j, int k = 0) const;
    // const T& at(const Vector3& pos) const;
    const T& at(const Vector3i& pos) const;
    const T& at(size_t i) const;
    const T& operator()(size_t x, size_t y, size_t z = 0) const;
    const T& operator()(size_t i) const;
    // const T& operator()(const Vector3& pos) const;
    const T& operator()(const Vector3i& pos) const;
    const T& operator[](size_t i) const;
    // const T& operator[](const Vector3& pos) const;
    const T& operator[](const Vector3i& pos) const;
    T& at(int i, int j, int k = 0);
    // T& at(const Vector3& pos);
    T& at(const Vector3i& pos);
    T& at(size_t i);
//    T& operator[](size_t x, size_t y);
    T& operator[](size_t i);
    // T& operator[](const Vector3& pos);
    T& operator[](const Vector3i& pos);
    T& operator()(size_t x, size_t y, size_t z = 0);
    T& operator()(size_t i);
    // T& operator()(const Vector3& pos);
    T& operator()(const Vector3i& pos);
    inline T& unsafe(size_t i) noexcept { return data[i]; }
    inline const T& unsafe(size_t i) const noexcept { return data[i]; }
    inline T& unsafe(size_t x, size_t y, size_t z = 0) noexcept { return data[getIndex(x, y, z)]; }
    inline const T& unsafe(size_t x, size_t y, size_t z = 0) const noexcept { return data[getIndex(x, y, z)]; }
    // inline T& unsafe(const Vector3& p) noexcept { return data[getIndex(p.x(), p.y(), p.z())]; }
    inline T& unsafe(const Vector3i& p) noexcept { return data[getIndex(p.x(), p.y(), p.z())]; }
    // inline const T& unsafe(const Vector3& p) const noexcept { return data[getIndex(p.x(), p.y(), p.z())]; }
    inline const T& unsafe(const Vector3i& p) const noexcept { return data[getIndex(p.x(), p.y(), p.z())]; }


    Vector3i getDimensions() const;
    int width() const;
    int depth() const;
    int height() const;

    int rows() const;
    int cols() const;

    inline size_t getIndex(size_t x, size_t y, size_t z) const noexcept;
    inline size_t getIndex(const Vector3i& coord) const noexcept;
    inline std::tuple<size_t, size_t, size_t> getCoord(size_t index) const noexcept;
    inline Vector3 getCoordAsVector3(size_t index) const;
    inline Vector3i getCoordAsVector3i(size_t index) const;
    inline bool checkCoord(int x, int y, int z = 0) const;
    inline bool checkCoord(const Vector3i& pos) const;
    inline bool checkIndex(size_t i) const;

    Matrix3<T> col(int colIndex, int depthIndex = 0);
    Matrix3<T> row(int rowIndex, int depthIndex = 0);

    template <class Func>
    void iterate(Func function) const;
    template <class Func>
    void iterateReverse(Func function) const;
    template <class Func>
    void iterateParallel(Func function) const;
    template <class Func>
    void iterateRandomly(Func function) const;

    T interpolate(const Vector3& coord, RETURN_VALUE_ON_OUTSIDE padding = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE) const;
    T interpolate(float x, float y, float z = 0, RETURN_VALUE_ON_OUTSIDE padding = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE) const;
    Matrix3<T>& addValueAt(T value, const Vector3& coord);
    Matrix3<T>& addValueAt(T value, float x, float y, float z = 0.f);

    int getNumberNeighbors(size_t x, size_t y, size_t z, bool using4connect = true) const;
    // int getNumberNeighbors(const Vector3& pos, bool using4connect = true) const;
    int getNumberNeighbors(const Vector3i& pos, bool using4connect = true) const;

    Matrix3<T> resize(float factor, RESIZE_MODE mode = LINEAR) const;
    Matrix3<T> resize(size_t newX, size_t newY, size_t newZ, RESIZE_MODE mode = LINEAR) const;
    // Matrix3<T> resize(const Vector3& newSize, RESIZE_MODE mode = LINEAR) const;
    Matrix3<T> resize(const Vector3i& newSize, RESIZE_MODE mode = LINEAR) const;

    Matrix3<T> resizeNearest(float factor) const;
    Matrix3<T> resizeNearest(size_t newX, size_t newY, size_t newZ) const;
    // Matrix3<T> resizeNearest(const Vector3& newSize) const;
    Matrix3<T> resizeNearest(const Vector3i& newSize) const;

    Matrix3<T> subset(int startX, int endX, int startY, int endY, int startZ = 0, int endZ = -1) const;
    // Matrix3<T> subset(const Vector3& start, const Vector3& end) const;
    Matrix3<T> subset(const Vector3i& start, const Vector3i& end) const;
    // Matrix3<T>& paste(const Matrix3<T>& matrixToPaste, const Vector3& upperLeftFrontCorner = Vector3());
    inline Matrix3<T>& paste(const Matrix3<T>& matrixToPaste, const Vector3i& upperLeftFrontCorner = Vector3i());
    Matrix3<T>& paste(const Matrix3<T>& matrixToPaste, int left, int up, int front);
    Matrix3<T>& add(const Matrix3<T>& matrixToAdd, const Vector3& upperLeftFrontCorner, bool useInterpolation = false);
    Matrix3<T>& add(const Matrix3<T>& matrixToAdd, const Vector3i& upperLeftFrontCorner);
    Matrix3<T>& add(const Matrix3<T>& matrixToAdd, int left, int up, int front, bool useInterpolation = false);
    Matrix3<T> concat(const Matrix3<T>& matrixToConcat);

    Matrix3<float> toDistanceMap(bool ignoreZlayer = false, bool considerBorders = true);
    Matrix3<std::complex<float>> FFT() const;
    Matrix3<std::complex<float>> iFFT(const Vector3i& finalDimensions = Vector3i::invalid) const;
    Matrix3<std::complex<float>> iFFT(size_t cropX, size_t cropY, size_t cropZ) const;

    Matrix3<T> flip(bool onX, bool onY = false, bool onZ = false);

    template <class U>
    Matrix3<T> convolution(const Matrix3<U>& convMatrix, CONVOLUTION_BORDERS border = ZERO_PAD) const;

    T min() const;
    T max() const;

    Matrix3<T>& min(const T& minVal);
    Matrix3<T>& max(const T& maxVal);

    // Matrix3<T>& max(const Matrix3<T>& otherMatrix, const Vector3& upperLeftFrontCorner);
    inline Matrix3<T>& max(const Matrix3<T>& otherMatrix, const Vector3i& upperLeftFrontCorner);
    Matrix3<T>& max(const Matrix3<T>& otherMatrix, int left, int up, int front);
    // Matrix3<T>& min(const Matrix3<T>& otherMatrix, const Vector3& upperLeftFrontCorner);
    inline Matrix3<T>& min(const Matrix3<T>& otherMatrix, const Vector3i& upperLeftFrontCorner);
    Matrix3<T>& min(const Matrix3<T>& otherMatrix, int left, int up, int front);

    static Matrix3<T> max(const Matrix3<T>& m1, const Matrix3<T>& m2);
    static Matrix3<T> min(const Matrix3<T>& m1, const Matrix3<T>& m2);

    Matrix3<T> abs() const;
    T sum() const;

    Matrix3<T> rounded(int precision = 0) const;

    Matrix3<T>& normalize();
    Matrix3<T> normalized() const;

    Matrix3<T>& normalizeUsing(NORMALIZE_METHOD normalizeMethod = NORMALIZE_MINMAX);
    Matrix3<T> normalizedUsing(NORMALIZE_METHOD normalizeMethod = NORMALIZE_MINMAX) const;

    Matrix3<T> transposeXY();

    Matrix3<Vector3> gradient() const;
    inline Matrix3<Vector3> grad() const;
    Matrix3<float> divergence() const;
    inline Matrix3<float> div() const;
    Matrix3<T> curl(float radius = 1.f) const;
    inline Matrix3<T> rot() const;
    Matrix3<T> laplacian() const;
    Vector3 gradient(const Vector3& position) const;
    inline Vector3 gradient(float posX, float posY, float posZ = 0) const;

    Matrix3<int> skeletonize() const;
    std::vector<CatmullRomSpline> skeletonizeToBSplines() const;
    Matrix3<Vector3> fillWithBSplines(std::vector<CatmullRomSpline> splines) const;
    Matrix3<T> dilate(bool use2D = false, float t = 1.f) const;
    Matrix3<T> erode(bool use2D = false, float t = 1.f) const;
    Matrix3<int> computeConnectedComponents(bool use4Connect = false) const;
    Matrix3<int> fillHoles(bool ignoreZlayer = true) const;
    Matrix3<int> findContour(bool use2D = false) const;
    std::vector<ShapeCurve> findContoursAsCurves() const; // 2D specific

    T trace() const;

    Matrix3<T> warpWith(const Matrix3<Vector3>& warper, int precision = 1) const;
    Matrix3<T> warpWith(const CatmullRomSpline& original, const CatmullRomSpline& warperCurve) const;

    Matrix3<T> warpWithoutInterpolation(const Matrix3<Vector3>& warper) const;
    Matrix3<T> warpWithoutInterpolation(const CatmullRomSpline& original, const CatmullRomSpline& warperCurve) const;

    static Matrix3<float> fbmNoise1D(FastNoiseLite noise, int sizeX, int sizeY, int sizeZ = 1);
    static Matrix3<Vector3> fbmNoise2D(FastNoiseLite noise, int sizeX, int sizeY, int sizeZ = 1);
    static Matrix3<Vector3> fbmNoise3D(FastNoiseLite noise, int sizeX, int sizeY, int sizeZ = 1);

    static Matrix3<float> gaussian(int sizeOnX, int sizeOnY, int sizeOnZ, float sigma, const Vector3& offset = Vector3());
    static Matrix3<float> normalizedGaussian(int sizeOnX, int sizeOnY, int sizeOnZ, float sigma, const Vector3& offset = Vector3());
    Matrix3<T> LaplacianOfGaussian(int sizeOnX, int sizeOnY, int sizeOnZ, float sigma) const;
    Matrix3<T> meanSmooth(int sizeOnX = 3, int sizeOnY = 3, int sizeOnZ = 3, [[maybe_unused]] bool ignoreBorders = false) const;
    Matrix3<T> gaussianSmooth(float sigma, bool ignoreZ = false, bool ignoreBorders = false, float limitFactor = 4.f) const;
    Matrix3<T> medianBlur(int sizeOnX = 3, int sizeOnY = 3, int sizeOnZ = 3, bool ignoreBorders = false) const;

    Matrix3<T>& insertRow(size_t indexToInsert, int affectedDimension, T newData = T());

    void clear() { this->sizeX = 0; this->sizeY = 0; this->sizeZ = 0; this->_updateStrides(); return this->data.clear(); }
    void reset(T newVal = T()) { for (auto& val : data) val = newVal; }

    Matrix3<int> binarize(T limitValue = T(), bool greaterValuesAreSetToOne = true, bool useAlsoTheEqualSign = false) const;
    Matrix3<int> binarizeBetween(T minValue, T maxValue, bool insideValuesAreSetToOne = true, bool useAlsoTheEqualSign = false) const;
    Matrix3<int> isosurface(T isovalue = T(), bool ignoreZtopBorder = false, bool ignoreBorders = true) const;

    Matrix3<T> slice(int index, int axis) const;
    Matrix3<T> sliceXY(int index) const;
    Matrix3<T> sliceYZ(int index) const;
    Matrix3<T> sliceXZ(int index) const;

    template <class U>
    operator Matrix3<U>() const {
        Matrix3<U> returned(this->getDimensions());
        for (size_t i = 0; i < this->size(); i++)
            returned[i] = (U)(this->data[i]);
        return returned;
    }


    // static Matrix3<T> random(const Vector3& dimensions);
    static Matrix3<T> random(const Vector3i& dimensions);
    static Matrix3<T> random(size_t sizeX, size_t sizeY, size_t sizeZ = 1);
    static Matrix3<T> identity(size_t sizeX, size_t sizeY, size_t sizeZ = 1);
    // static Matrix3<T> perlin(const Vector3& dimensions, const Vector3 &scale = Vector3(1, 1, 1), int seed = 0);
    static Matrix3<T> perlin(const Vector3i& dimensions, const Vector3 &scale = Vector3(1, 1, 1), int seed = 0);

    static Matrix3<Vector3> fromImageRGB(const std::string& filename);
    static Matrix3<float> fromImageBW(const std::string& filename);

    Matrix3<T> operator-() const;
    template <class U>
    Matrix3<T>& operator+=(const Matrix3<U>& o);
    template <class U>
    Matrix3<T>& operator-=(const Matrix3<U>& o);
    template <class U>
    Matrix3<T>& operator*=(const Matrix3<U>& o);
    template <class U>
    Matrix3<T>& operator/=(const Matrix3<U>& o);
    template <class U>
    Matrix3<T>& operator*=(U o);
    template <class U>
    Matrix3<T>& operator/=(U o);
    template <class U>
    Matrix3<T>& operator+=(U o);
    template <class U>
    Matrix3<T>& operator-=(U o);
    template <class U>
    bool operator==(const Matrix3<U>& o) const;


    std::string toString() const {return "Matrix3 (" + std::to_string(this->sizeX) + "x" + std::to_string(this->sizeY) + "x" + std::to_string(this->sizeZ) + ")"; }

    std::vector<T> data;
    std::size_t sizeX = 0, sizeY = 0, sizeZ = 1;
    T dummyValue = T(); // Trash value, just for invalid at() calls

    bool raiseErrorOnBadCoord = false; // THIS SHOULD CLEARLY BE SET TO TRUE, BUT F**K IT!
    T defaultValueOnBadCoord = T();
    RETURN_VALUE_ON_OUTSIDE returned_value_on_outside = DEFAULT_VALUE;
    bool stillRaiseErrorForX = false;
    bool stillRaiseErrorForY = false;
    bool stillRaiseErrorForZ = false;

    inline auto begin() const noexcept{ return data.begin(); }
    inline auto end() const noexcept{ return data.end(); }
    inline auto begin() { return data.begin(); }
    inline auto end() { return data.end(); }
    inline std::size_t size() const noexcept{ return end() - begin(); }
    inline bool empty() const noexcept{ return begin() == end(); }

    inline Vector3 getMirrorPosition(const Vector3& pos) const noexcept;
    inline Vector3i getMirrorPosition(const Vector3i& pos) const noexcept;
    inline size_t getMirrorIndex(int x, int y, int z = 0) const noexcept;
    inline size_t getMirrorIndex(const Vector3i& p) const noexcept { return getMirrorIndex(p.x(), p.y(), p.z()); }
    inline Vector3 getWrappedPosition(const Vector3& pos) const noexcept;
    inline Vector3i getWrappedPosition(const Vector3i& pos) const noexcept;
    inline size_t getWrappedIndex(int x, int y, int z = 0) const noexcept;
    inline size_t getWrappedIndex(const Vector3i& p) const noexcept { return getWrappedIndex(p.x(), p.y(), p.z()); }
    inline Vector3 getRepeatPosition(const Vector3& pos) const noexcept;
    inline Vector3i getRepeatPosition(const Vector3i& pos) const noexcept;
    inline size_t getRepeatIndex(int x, int y, int z = 0) const noexcept;
    inline size_t getRepeatIndex(const Vector3i& p) const noexcept { return getRepeatIndex(p.x(), p.y(), p.z()); }

    Matrix3& init(const std::vector<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ);

    std::string displayValues() const;
    std::string displayAsPlot(T min = 0.f, T max = 0.f, std::vector<std::string> patterns = {}, std::map<T, std::string> specialCharactersAtValue = {}, T specialCharEpsilon = 1e-5, std::string charForError = "X", std::string separator = "") const;

    void _updateStrides() {
        _strideY = sizeX;
        _strideZ = _strideY * sizeY;
    }
// protected:
    std::size_t _strideY = 0, _strideZ = 0;
};

template <class T> template <class Func>
void Matrix3<T>::iterate(Func function) const
{
    for (size_t z = 0; z < this->sizeZ; z++) {
        for (size_t y = 0; y < this->sizeY; y++) {
                for (size_t x = 0; x < this->sizeX; x++) {
                if constexpr (std::is_invocable_v<Func, Vector3i>) {
                    Vector3i pos(x, y, z);
                    function(pos);
                } else if constexpr (std::is_invocable_v<Func, int, int, int>) {
                    Vector3i pos(x, y, z);
                    function(pos.x(), pos.y(), pos.z());
                } else if constexpr (std::is_invocable_v<Func, int>) {
                    function(this->getIndex(x, y, z));
                } else {
                    function();
                }
            }
        }
    }
    /*
    for (size_t i = 0; i < this->size(); i++) {
        if constexpr (std::is_invocable_v<Func, Vector3>) {
            Vector3 pos = this->getCoordAsVector3(i);
            function(pos);
        } else if constexpr (std::is_invocable_v<Func, int, int, int>) {
            Vector3 pos = this->getCoordAsVector3(i);
            function(pos.x, pos.y, pos.z);
        } else if constexpr (std::is_invocable_v<Func, int>) {
            function(i);
        } else {
            function();
        }
    }*/
}
template <class T> template <class Func>
void Matrix3<T>::iterateReverse(Func function) const
{
    for (int z = this->sizeZ - 1; z >= 0; z--) {
        for (int y = this->sizeY - 1; y >= 0; y--) {
                for (int x = this->sizeX - 1; x >= 0; x--) {
                if constexpr (std::is_invocable_v<Func, Vector3i>) {
                    Vector3i pos(x, y, z);
                    function(pos);
                } else if constexpr (std::is_invocable_v<Func, int, int, int>) {
                    Vector3i pos(x, y, z);
                    function(pos.x(), pos.y(), pos.z());
                } else if constexpr (std::is_invocable_v<Func, int>) {
                    function(this->getIndex(x, y, z));
                } else {
                    function();
                }
            }
        }
    }
    /*
    for (int i = int(this->size()) - 1; i >= 0; i--) {
        if constexpr (std::is_invocable_v<Func, Vector3>) {
            Vector3 pos = this->getCoordAsVector3(i);
            function(pos);
        } else if constexpr (std::is_invocable_v<Func, int, int, int>) {
            Vector3 pos = this->getCoordAsVector3(i);
            function(pos.x, pos.y, pos.z);
        } else if constexpr (std::is_invocable_v<Func, int>) {
            function(i);
        } else {
            function();
        }
    }
    */
}
template <class T> template <class Func>
void Matrix3<T>::iterateParallel(Func function) const
{
    if constexpr (std::is_invocable_v<Func, int>) {
        #pragma omp parallel for
        for (size_t i = 0; i < this->size(); i++) {
            function(i);
        }
    } else {
        #pragma omp parallel for collapse(3) schedule(dynamic)
        for (size_t z = 0; z < this->sizeZ; z++) {
            for (size_t y = 0; y < this->sizeY; y++) {
                for (size_t x = 0; x < this->sizeX; x++) {
                    if constexpr (std::is_invocable_v<Func, Vector3i>) {
                        Vector3i pos(x, y, z);
                        function(pos);
                    } else if constexpr (std::is_invocable_v<Func, int, int, int>) {
                        Vector3i pos(x, y, z);
                        function(pos.x(), pos.y(), pos.z());
                    } else /*if constexpr (std::is_invocable_v<Func, int>) {
                        function(this->getIndex(x, y, z));
                    } else*/ {
                        function();
                    }
                }
            }
        }
    }
    /*
    #pragma omp parallel for
    for (size_t i = 0; i < this->size(); i++) {
        if constexpr (std::is_invocable_v<Func, Vector3>) {
            Vector3 pos = this->getCoordAsVector3(i);
            function(pos);
        } else if constexpr (std::is_invocable_v<Func, int, int, int>) {
            Vector3 pos = this->getCoordAsVector3(i);
            function(pos.x, pos.y, pos.z);
        } else if constexpr (std::is_invocable_v<Func, int>) {
            function(i);
        } else {
            function();
        }
    }
    */
}
template <class T> template <class Func>
void Matrix3<T>::iterateRandomly(Func function) const
{
    std::vector<size_t> iter(this->size());
    for (size_t i = 0; i < iter.size(); i++)
        iter[i] = i;
    std::shuffle(iter.begin(), iter.end(), random_gen::random_generator);
    for (size_t j = 0; j < iter.size(); j++) {
        size_t i = iter[j];
        if constexpr (std::is_invocable_v<Func, Vector3i>) {
            Vector3i pos = this->getCoordAsVector3i(i);
            function(pos);
        } else if constexpr (std::is_invocable_v<Func, int, int, int>) {
            Vector3i pos = this->getCoordAsVector3i(i);
            function(pos.x(), pos.y(), pos.z());
        } else if constexpr (std::is_invocable_v<Func, int>) {
            function(i);
        } else {
            function();
        }
    }
}
//template <class T>
Matrix3<float> operator-(const float a, const Matrix3<float>& b);
//template <class T>
Matrix3<float> operator+(const float a, const Matrix3<float>& b);
#include <sstream>
template <class T>
std::string Matrix3<T>::displayValues() const
{
    std::stringstream out;
    for (size_t z = 0; z < this->sizeZ; z++) {
        out << "[Z-level = " << z << "] : \n";
        for (size_t y = 0; y < this->sizeY; y++) {
            for (size_t x = 0; x < this->sizeX; x++) {
                out << std::setw(2) << at(x, y, z) << "\t";
            }
            out << "\n";
        }
    }
    return out.str();
}

template <class T>
std::string Matrix3<T>::displayAsPlot(T min, T max, std::vector<std::string> patterns, std::map<T, std::string> specialCharactersAtValue, T specialCharEpsilon, std::string charForError, std::string separator) const
{
    if (patterns.empty())
        patterns = {".", "-", "=", "#"};

    if (min == 0.f && max == 0.f) {
        min = this->min();
        max = this->max();
    }
    std::stringstream out;
    for (size_t z = 0; z < this->sizeZ; z++) {
        out << "[Z-level = " << z << "] : \n";
        for (size_t y = 0; y < this->sizeY; y++) {
            for (size_t x = 0; x < this->sizeX; x++) {
                T val = this->unsafe(x, y, z);
                bool specialCharUsed = false;
                for (const auto& [targetValue, character] : specialCharactersAtValue) {
                    if (targetValue - specialCharEpsilon <= val && val <= targetValue + specialCharEpsilon) {
                        out << character << separator;
                        specialCharUsed = true;
                        break; // Don't add other characters
                    }
                }
                if (!specialCharUsed) {
                    float prop = interpolation::linear(val, min, max);
                    if (prop < 0.f || prop > 1.f || isnan(prop))
                        out << charForError << separator;
                    else
                        out << patterns[int(prop * float(patterns.size() - 1))] << separator;
                }
            }
            out << "\n";
        }
    }
    return out.str();
}



template <class T>
Matrix3<T>::Matrix3()
{
}
/*template <class T>
Matrix3<T>::Matrix3(const Vector3& size, T initValue) : Matrix3<T>(size.x(), size.y(), size.z(), initValue)
{
}*/
template <class T>
Matrix3<T>::Matrix3(const Vector3i& size, T initValue) : Matrix3<T>(size.x(), size.y(), size.z(), initValue)
{
}
template <class T>
Matrix3<T>::Matrix3(size_t sizeX, size_t sizeY, size_t sizeZ, T initValue)
{
    this->data = std::vector<T>(sizeX * sizeY * sizeZ, initValue);
    this->sizeX = sizeX;
    this->sizeY = sizeY;
    this->sizeZ = sizeZ;
    _updateStrides();
}
template <class T>
Matrix3<T>::Matrix3(const std::vector<T>& data, size_t sizeX, size_t sizeY, int sizeZ)
{
    if (sizeZ == -1) {
        sizeZ = int(data.size()) / (sizeX * sizeY);
    }
    init(data, sizeX, sizeY, sizeZ);
}
template <class T>
Matrix3<T>::Matrix3(const std::vector<std::vector<T>>& data)
{
    std::vector<T> oneMatrix;
    for (const std::vector<T>& row : data)
        oneMatrix.insert(oneMatrix.end(), row.begin(), row.end());
    int sizeX = data[0].size();
    int sizeY = data.size();
    init(oneMatrix, sizeX, sizeY, 1);
}
template <class T>
Matrix3<T>::Matrix3(const std::vector<std::vector<std::vector<T> > > &data)
{
    std::vector<T> oneMatrix;
    for (std::vector<std::vector<T>>& grid : data)
        for (std::vector<T>& row : grid)
            oneMatrix.insert(oneMatrix.end(), row.begin(), row.end());
    int sizeX = data[0].size();
    int sizeY = data.size();
    int sizeZ = (sizeY > 0 ? data[0][0].size() : 0);
    init(oneMatrix, sizeX, sizeY, sizeZ);
}

template <class T>
inline bool Matrix3<T>::checkCoord(int x, int y, int z) const
{
    return ((0 <= x && static_cast<std::size_t>(x) < sizeX) && (0 <= y && static_cast<std::size_t>(y) < sizeY) && (0 <= z && static_cast<std::size_t>(z) < sizeZ));
}

template <class T>
inline bool Matrix3<T>::checkCoord(const Vector3i& pos) const
{
    // if (pos.minComp() < 0 || pos.x() > sizeX-1 || pos.y() > sizeY-1 || pos.z() > sizeZ-1) return false;
    return checkCoord(pos.x(), pos.y(), pos.z());
}

template <class T>
inline bool Matrix3<T>::checkIndex(size_t i) const
{
    return (0 <= i && i < data.size());
}

template <class T>
Matrix3<T> Matrix3<T>::col(int colIndex, int depthIndex)
{
    Matrix3<T> res(1, this->sizeY, 1);
    auto dst = res.data.data();
    const auto src = data.data();

    #pragma omp parallel for
    for (size_t y = 0; y < this->sizeY; y++) {
        dst[y] = src[getIndex(colIndex, y, depthIndex)];
    }
    return res;
}

template <class T>
Matrix3<T> Matrix3<T>::row(int rowIndex, int depthIndex)
{
    Matrix3<T> res(this->sizeX, 1, 1);
    auto dst = res.data.data();
    const auto src = data.data();
    #pragma omp parallel for
    for (size_t x = 0; x < this->sizeX; x++) {
        dst[x] = src[getIndex(x, rowIndex, depthIndex)];
    }
    return res;
}

template <class T>
T Matrix3<T>::interpolate(const Vector3& coord, RETURN_VALUE_ON_OUTSIDE padding) const
{
    Vector3i round = coord.floor();
    Vector3 cellOffset = coord - round;

    using PosFn = size_t (Matrix3::*)(const Vector3i&) const;
    PosFn pos_fn = nullptr;
    switch (padding) {
    case RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE:
        pos_fn = &Matrix3::getMirrorIndex; break;
    case RETURN_VALUE_ON_OUTSIDE::WRAPPED_VALUE:
        pos_fn = &Matrix3::getWrappedIndex; break;
    case RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE:
        pos_fn = &Matrix3::getRepeatIndex; break;
    default:
        break;
    }

    auto mapPos = [&](const Vector3i& p) -> size_t {
        return pos_fn ? (this->*pos_fn)(p) : getIndex(p);
    };

    const auto src = data.data();

    // bool previousErrorConfig = this->raiseErrorOnBadCoord;
    // RETURN_VALUE_ON_OUTSIDE previousOutsideConfig = this->returned_value_on_outside;
//    this->raiseErrorOnBadCoord = false;
//    this->returned_value_on_outside = padding;
    T f000 = src[mapPos(round + Vector3i(0, 0, 0))];
    T f100 = src[mapPos(round + Vector3i(1, 0, 0))];
    T f010 = src[mapPos(round + Vector3i(0, 1, 0))];
    T f110 = src[mapPos(round + Vector3i(1, 1, 0))];
    T f001 = src[mapPos(round + Vector3i(0, 0, 1))];
    T f101 = src[mapPos(round + Vector3i(1, 0, 1))];
    T f011 = src[mapPos(round + Vector3i(0, 1, 1))];
    T f111 = src[mapPos(round + Vector3i(1, 1, 1))];
    // Interpolation
    T interpol = ((
                              f000 * (1-cellOffset.x()) + f100 * cellOffset.x()) * (1-cellOffset.y()) + (
                              f010 * (1-cellOffset.x()) + f110 * cellOffset.x()) * cellOffset.y()) * (1 - cellOffset.z()) +
                        ((
                             f001 * (1-cellOffset.x()) + f101 * cellOffset.x()) * (1-cellOffset.y()) + (
                             f011 * (1-cellOffset.x()) + f111 * cellOffset.x()) * cellOffset.y()) * cellOffset.z();
//    this->raiseErrorOnBadCoord = previousErrorConfig;
//    this->returned_value_on_outside = previousOutsideConfig;
    return interpol;
}

template <class T>
T Matrix3<T>::interpolate(float x, float y, float z, RETURN_VALUE_ON_OUTSIDE padding) const
{
    return interpolate(Vector3(x, y, z), padding);
}

/*template <class T>
const T &Matrix3<T>::at(const Vector3& pos) const
{
    return this->at(pos.x(), pos.y(), pos.z());
}*/
template <class T>
const T &Matrix3<T>::at(const Vector3i& pos) const
{
    return this->at(pos.x(), pos.y(), pos.z());
}
template <class T>
const T &Matrix3<T>::at(int i, int j, int k) const
{
    if (checkCoord(i, j, k)) {
        int index = getIndex(i, j, k);
        return this->data[index];
    }
    bool raiseError = raiseErrorOnBadCoord;
    Vector3i newPos(i, j, k);
    if (!raiseErrorOnBadCoord) {
        if (returned_value_on_outside == DEFAULT_VALUE)
            return defaultValueOnBadCoord;

        if (stillRaiseErrorForX && (i < 0 || sizeX <= static_cast<std::size_t>(i)))
            return defaultValueOnBadCoord;
        if (stillRaiseErrorForY && (j < 0 || sizeY <= static_cast<std::size_t>(j)))
            return defaultValueOnBadCoord;
        if (stillRaiseErrorForZ && (k < 0 || sizeZ <= static_cast<std::size_t>(k)))
            return defaultValueOnBadCoord;

        if (returned_value_on_outside == MIRROR_VALUE)
            newPos = getMirrorPosition(newPos);
        else if (returned_value_on_outside == WRAPPED_VALUE)
            newPos = getWrappedPosition(newPos);
        else if (returned_value_on_outside == REPEAT_VALUE)
            newPos = getRepeatPosition(newPos);

    }
    if (!raiseError) {
        return this->unsafe(newPos.x(), newPos.y(), newPos.z());
    }
    else
        throw std::out_of_range("Trying to access coord (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ") on matrix of size "
            + std::to_string(sizeX) + "x" + std::to_string(sizeY) + "x" + std::to_string(sizeZ) + ". Max index is " + std::to_string(sizeX * sizeY * sizeZ - 1));
}

template <class T>
const T &Matrix3<T>::at(size_t i) const
{
    if (i < data.size()) {
        return this->data[i];
    }
    int x, y, z;
    std::tie(x, y, z) = this->getCoord(i);

    bool raiseError = raiseErrorOnBadCoord;
    Vector3i newPos(x, y, z);
    if (!raiseErrorOnBadCoord) {
        if (returned_value_on_outside == DEFAULT_VALUE)
            return defaultValueOnBadCoord;

        if (stillRaiseErrorForX && (x < 0 || sizeX <= static_cast<std::size_t>(x)))
            return defaultValueOnBadCoord;
        if (stillRaiseErrorForY && (y < 0 || sizeY <= static_cast<std::size_t>(y)))
            return defaultValueOnBadCoord;
        if (stillRaiseErrorForZ && (z < 0 || sizeZ <= static_cast<std::size_t>(z)))
            return defaultValueOnBadCoord;

        if (returned_value_on_outside == MIRROR_VALUE)
            newPos = getMirrorPosition(newPos);
        else if (returned_value_on_outside == WRAPPED_VALUE)
            newPos = getWrappedPosition(newPos);
        else if (returned_value_on_outside == REPEAT_VALUE)
            newPos = getRepeatPosition(newPos);
    }
    if (!raiseError)
        return this->unsafe(newPos.x(), newPos.y(), newPos.z());
    else
        throw std::out_of_range("Trying to access index " + std::to_string(i) + " (coord " + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ") on matrix of size "
            + std::to_string(sizeX) + "x" + std::to_string(sizeY) + "x" + std::to_string(sizeZ) + ". Max index is " + std::to_string(sizeX * sizeY * sizeZ - 1));
}
template <class T>
const T& Matrix3<T>::operator()(size_t x, size_t y, size_t z) const {
    return this->at(x, y, z);
}
template <class T>
const T& Matrix3<T>::operator()(size_t i) const {
    return this->at(i);
}
/*template <class T>
const T& Matrix3<T>::operator()(const Vector3& pos) const {
    return this->at(pos);
}*/
template <class T>
const T& Matrix3<T>::operator()(const Vector3i& pos) const {
    return this->at(pos);
}
template <class T>
const T& Matrix3<T>::operator[](size_t i) const {
    return this->at(i);
}
/*template <class T>
const T& Matrix3<T>::operator[](const Vector3& pos) const {
    return this->at(pos);
}*/
template <class T>
const T& Matrix3<T>::operator[](const Vector3i& pos) const {
    return this->at(pos);
}

/*template <class T>
T &Matrix3<T>::at(const Vector3& pos)
{
    return this->at(pos.x(), pos.y(), pos.z());
}*/
template <class T>
T &Matrix3<T>::at(const Vector3i& pos)
{
    return this->at(pos.x(), pos.y(), pos.z());
}
template <class T>
T &Matrix3<T>::at(int i, int j, int k)
{
    if (this->empty()) {
        throw std::out_of_range("Grid is empty : " + this->toString());
    }
    if (checkCoord(i, j, k)) {
        int index = getIndex(i, j, k);
        return this->data[index];
    }
    bool raiseError = raiseErrorOnBadCoord;
    dummyValue = defaultValueOnBadCoord;
    Vector3i newPos(i, j, k);
    if (!raiseErrorOnBadCoord) {
        if (returned_value_on_outside == DEFAULT_VALUE)
            return dummyValue; // defaultValueOnBadCoord;

        if (stillRaiseErrorForX && (i < 0 || sizeX <= static_cast<std::size_t>(i)))
            return dummyValue; // defaultValueOnBadCoord;
        if (stillRaiseErrorForY && (j < 0 || sizeY <= static_cast<std::size_t>(j)))
            return dummyValue; // defaultValueOnBadCoord;
        if (stillRaiseErrorForZ && (k < 0 || sizeZ <= static_cast<std::size_t>(k)))
            return dummyValue; // defaultValueOnBadCoord;

        if (returned_value_on_outside == MIRROR_VALUE)
            newPos = getMirrorPosition(newPos);
        else if (returned_value_on_outside == WRAPPED_VALUE)
            newPos = getWrappedPosition(newPos);
        else if (returned_value_on_outside == REPEAT_VALUE)
            newPos = getRepeatPosition(newPos);

    }
    if (!raiseError)
        return this->unsafe(newPos.x(), newPos.y(), newPos.z());
    else
        throw std::out_of_range("Trying to access coord (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ") on matrix of size "
            + std::to_string(sizeX) + "x" + std::to_string(sizeY) + "x" + std::to_string(sizeZ) + ". Max index is " + std::to_string(sizeX * sizeY * sizeZ - 1));
}

template <class T>
T &Matrix3<T>::at(size_t i)
{
    if (i >= 0 && i < sizeX * sizeY * sizeZ) {
        return this->data[i];
    }
    int x, y, z;
    std::tie(x, y, z) = this->getCoord(i);

    bool raiseError = raiseErrorOnBadCoord;
    Vector3i newPos(x, y, z);
    if (!raiseErrorOnBadCoord) {
        if (returned_value_on_outside == DEFAULT_VALUE)
            return dummyValue; // defaultValueOnBadCoord;

        if (stillRaiseErrorForX && (x < 0 || sizeX <= static_cast<std::size_t>(x)))
            return dummyValue; // defaultValueOnBadCoord;
        if (stillRaiseErrorForY && (y < 0 || sizeY <= static_cast<std::size_t>(y)))
            return dummyValue; // defaultValueOnBadCoord;
        if (stillRaiseErrorForZ && (z < 0 || sizeZ <= static_cast<std::size_t>(z)))
            return dummyValue; // defaultValueOnBadCoord;

        if (returned_value_on_outside == MIRROR_VALUE)
            newPos = getMirrorPosition(newPos);
        else if (returned_value_on_outside == WRAPPED_VALUE)
            newPos = getWrappedPosition(newPos);
        else if (returned_value_on_outside == REPEAT_VALUE)
            newPos = getRepeatPosition(newPos);
    }
    if (!raiseError)
        return this->unsafe(newPos.x(), newPos.y(), newPos.z());
    else
        throw std::out_of_range("Trying to access index " + std::to_string(i) + " (coord " + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ") on matrix of size "
            + std::to_string(sizeX) + "x" + std::to_string(sizeY) + "x" + std::to_string(sizeZ) + ". Max index is " + std::to_string(sizeX * sizeY * sizeZ - 1));
}
//template <class T>
//T& Matrix3<T>::operator[](size_t x, size_t y) {
//    return this->at(x, y);
//}
template <class T>
T& Matrix3<T>::operator[](size_t i) {
    return this->at(i);
}
/*template <class T>
T& Matrix3<T>::operator[](const Vector3& pos) {
    return this->at(pos);
}*/
template <class T>
T& Matrix3<T>::operator[](const Vector3i& pos) {
    return this->at(pos);
}
template <class T>
T& Matrix3<T>::operator()(size_t x, size_t y, size_t z) {
    return this->at(x, y, z);
}
template <class T>
T& Matrix3<T>::operator()(size_t i) {
    return this->at(i);
}
/*template <class T>
T& Matrix3<T>::operator()(const Vector3& pos) {
    return this->at(pos);
}*/
template <class T>
T& Matrix3<T>::operator()(const Vector3i& pos) {
    return this->at(pos);
}

template <class T>
inline size_t Matrix3<T>::getIndex(size_t x, size_t y, size_t z) const noexcept
{
    return z * this->_strideZ + y * this->_strideY + x;
}
template <class T>
inline size_t Matrix3<T>::getIndex(const Vector3i& coord) const noexcept
{
    return this->getIndex(coord.x(), coord.y(), coord.z());
}

template <class T>
inline std::tuple<size_t, size_t, size_t> Matrix3<T>::getCoord(size_t index) const noexcept
{
    const size_t z = index / _strideZ;
    index -= z * _strideZ;
    const size_t y = index / _strideY;
    const size_t x = index - y * _strideY;
    return {x, y, z};
}
template <class T>
inline Vector3 Matrix3<T>::getCoordAsVector3(size_t index) const
{
    const size_t z = index / _strideZ;
    index -= z * _strideZ;
    const size_t y = index / _strideY;
    const size_t x = index - y * _strideY;
    return Vector3(x, y, z);
}
template <class T>
inline Vector3i Matrix3<T>::getCoordAsVector3i(size_t index) const
{
    const size_t z = index / _strideZ;
    index -= z * _strideZ;
    const size_t y = index / _strideY;
    const size_t x = index - y * _strideY;
    return Vector3i(x, y, z);
}

template <class T>
Matrix3<T>& Matrix3<T>::addValueAt(T value, const Vector3& coord) {
    bool previousError = this->raiseErrorOnBadCoord;

    this->raiseErrorOnBadCoord = false;

    Vector3i floorPos = coord.floor();
    Vector3 offset = coord - floorPos;

    const Vector3i v0 = floorPos + Vector3i(0, 0, 0);
    const Vector3i v1 = floorPos + Vector3i(0, 0, 1);
    const Vector3i v2 = floorPos + Vector3i(0, 1, 0);
    const Vector3i v3 = floorPos + Vector3i(0, 1, 1);
    const Vector3i v4 = floorPos + Vector3i(1, 0, 0);
    const Vector3i v5 = floorPos + Vector3i(1, 0, 1);
    const Vector3i v6 = floorPos + Vector3i(1, 1, 0);
    const Vector3i v7 = floorPos + Vector3i(1, 1, 1);

    this->at(v0) += value * (1 - offset.x()) * (1 - offset.y()) * (1 - offset.z());
    this->at(v1) += value * (1 - offset.x()) * (1 - offset.y()) * (    offset.z());
    this->at(v2) += value * (1 - offset.x()) * (    offset.y()) * (1 - offset.z());
    this->at(v3) += value * (1 - offset.x()) * (    offset.y()) * (    offset.z());
    this->at(v4) += value * (    offset.x()) * (1 - offset.y()) * (1 - offset.z());
    this->at(v5) += value * (    offset.x()) * (1 - offset.y()) * (    offset.z());
    this->at(v6) += value * (    offset.x()) * (    offset.y()) * (1 - offset.z());
    this->at(v7) += value * (    offset.x()) * (    offset.y()) * (    offset.z());

    this->raiseErrorOnBadCoord = previousError;
    return *this;
}
template <class T>
Matrix3<T>& Matrix3<T>::addValueAt(T value, float x, float y, float z) {
    return this->addValueAt(value, Vector3(x, y, z));
}

template <class T>
Vector3 Matrix3<T>::gradient(const Vector3& position) const
{
    Vector3i p = position.floor();

    auto mirror = [&](const Vector3i& q) -> Vector3i {
        return this->checkCoord(q) ? q : this->getMirrorPosition(q);
    };

    const Vector3i p0 = mirror(p);
    const Vector3i px = mirror(p + Vector3i(1, 0, 0));
    const Vector3i py = mirror(p + Vector3i(0, 1, 0));
    const Vector3i pz = mirror(p + Vector3i(0, 0, 1));

    const T f0 = this->unsafe(p0);
    return Vector3(
        this->unsafe(px) - f0,
        this->unsafe(py) - f0,
        this->unsafe(pz) - f0
        );
    /*
    auto self = *this;
    self.raiseErrorOnBadCoord = false;
    self.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE;
    Vector3 flooredPos = position.floor();
//    Vector3 offset = position - flooredPos;

    return Vector3(
                (self.at(flooredPos + Vector3(1, 0, 0)) - self.at(flooredPos)), // * offset.x,
                (self.at(flooredPos + Vector3(0, 1, 0)) - self.at(flooredPos)), // * offset.y,
                (self.at(flooredPos + Vector3(0, 0, 1)) - self.at(flooredPos)) // * offset.z
                );
    */
}

template <class T>
Vector3 Matrix3<T>::gradient(float posX, float posY, float posZ) const
{
    return gradient(Vector3(posX, posY, posZ));
}

template <class T>
Matrix3<T> Matrix3<T>::dilate(bool use2D, float t) const
{
    if (this->empty()) return *this;
    if (t <= 0.f) return *this;

    const size_t sx = this->sizeX, sy = this->sizeY, sz = this->sizeZ;
    const size_t xy = sx * sy;

    // ping-pong buffers
    Matrix3<T> cur = *this;
    Matrix3<T> nxt(sx, sy, sz);

    auto clampi = [](int v, size_t lo, size_t hi) -> size_t {
        return v < int(lo) ? lo : (v > int(hi) ? hi : static_cast<std::size_t>(v));
    };

    // number of full iterations + one fractional blend at the end
    const int steps = (int)std::floor(t);
    const float frac = t - (float)steps;

    auto one_step = [&](float blend /* 1.0 for full step, frac for fractional */) {
        // Full morphological dilation gives maxVal.
        // If blend < 1: interpolate between original and dilated value
        #pragma omp parallel for collapse(3)
        for (size_t z = 0; z < sz; ++z) {
            for (size_t y = 0; y < sy; ++y) {
                for (size_t x = 0; x < sx; ++x) {

                    const size_t z0 = use2D ? z : clampi(z - 1, 0, sz - 1);
                    const size_t z1 = use2D ? z : clampi(z + 1, 0, sz - 1);

                    T maxVal = cur.data[z * xy + y * sx + x];

                    // iterate neighborhood with clamp borders (REPEAT_VALUE behavior)
                    for (size_t zz = z0; zz <= z1; ++zz) {
                        const size_t zbase = zz * xy;
                        for (size_t yy = clampi(y - 1, 0, sy - 1); yy <= clampi(y + 1, 0, sy - 1); ++yy) {
                            const size_t ybase = zbase + yy * sx;
                            const size_t x0 = clampi(x - 1, 0, sx - 1);
                            const size_t x1 = clampi(x + 1, 0, sx - 1);
                            for (size_t xx = x0; xx <= x1; ++xx) {
                                maxVal = std::max(maxVal, cur.data[ybase + xx]);
                            }
                        }
                    }

                    const T center = cur.data[z * xy + y * sx + x];
                    nxt.data[z * xy + y * sx + x] = (blend >= 1.f)
                                                        ? maxVal
                                                        : (T)((1.f - blend) * center + blend * maxVal);
                }
            }
        }

        cur.data.swap(nxt.data);
    };

    // full integer steps
    for (int i = 0; i < steps; ++i) one_step(1.f);
    // fractional step if needed
    if (frac > 0.f) one_step(frac);

    return cur;
}

template <class T>
Matrix3<T> Matrix3<T>::erode(bool use2D, float t) const
{
    if (this->empty()) return *this;
    if (t <= 0.f) return *this;

    const size_t sx = this->sizeX, sy = this->sizeY, sz = this->sizeZ;
    const size_t xy = sx * sy;

    // ping-pong buffers
    Matrix3<T> cur = *this;
    Matrix3<T> nxt(sx, sy, sz);

    auto clampi = [](int v, size_t lo, size_t hi) -> size_t {
        return v < int(lo) ? lo : (v > int(hi) ? hi : static_cast<std::size_t>(v));
    };

    // number of full iterations + one fractional blend at the end
    const int steps = (int)std::floor(t);
    const float frac = t - (float)steps;

    auto one_step = [&](float blend /* 1.0 for full step, frac for fractional */) {
        // Full morphological erosion gives minVal.
        // If blend < 1: interpolate between original and eroded value
        #pragma omp parallel for collapse(3)
        for (size_t z = 0; z < sz; ++z) {
            for (size_t y = 0; y < sy; ++y) {
                for (size_t x = 0; x < sx; ++x) {

                    const size_t z0 = use2D ? z : clampi(z - 1, 0, sz - 1);
                    const size_t z1 = use2D ? z : clampi(z + 1, 0, sz - 1);

                    T minVal = cur.data[z * xy + y * sx + x];

                    // iterate neighborhood with clamp borders (REPEAT_VALUE behavior)
                    for (size_t zz = z0; zz <= z1; ++zz) {
                        const size_t zbase = zz * xy;
                        for (size_t yy = clampi(y - 1, 0, sy - 1); yy <= clampi(y + 1, 0, sy - 1); ++yy) {
                            const size_t ybase = zbase + yy * sx;
                            const size_t x0 = clampi(x - 1, 0, sx - 1);
                            const size_t x1 = clampi(x + 1, 0, sx - 1);
                            for (size_t xx = x0; xx <= x1; ++xx) {
                                minVal = std::min(minVal, cur.data[ybase + xx]);
                            }
                        }
                    }

                    const T center = cur.data[z * xy + y * sx + x];
                    nxt.data[z * xy + y * sx + x] = (blend >= 1.f)
                                                        ? minVal
                                                        : (T)((1.f - blend) * center + blend * minVal);
                }
            }
        }

        cur.data.swap(nxt.data);
    };

    // full integer steps
    for (int i = 0; i < steps; ++i) one_step(1.f);
    // fractional step if needed
    if (frac > 0.f) one_step(frac);

    return cur;
}


template <class T>
Matrix3<int> Matrix3<T>::fillHoles(bool ignoreZlayer) const
{
    if (ignoreZlayer) {
        Matrix3<int> components(getDimensions() + Vector3(2, 2), 1);
        components.paste((1.f - *this), Vector3(1, 1));
        for (auto& c : components)
            c = clamp(c, 0, 1);
        components = ((Matrix3<int>)(components)).computeConnectedComponents(true);
        auto outsideValue = components[0];

        components.iterateParallel([&] (size_t i) {
            components[i] = (components[i] == outsideValue ? 0 : 1);
        });
        return components.subset(1, components.sizeX - 1, 1, components.sizeY - 1);
    } else {
        Matrix3<int> components(getDimensions() + Vector3(2, 2, 2), 1);
        components.paste((1.f - *this), Vector3(1, 1, 1));
        for (auto& c : components)
            c = clamp(c, 0, 1);
        components = ((Matrix3<int>)(components)).computeConnectedComponents(true);
        auto outsideValue = components[0];

        components.iterateParallel([&] (size_t i) {
            components[i] = (components[i] == outsideValue ? 0 : 1);
        });
        return components.subset(1, components.sizeX - 1, 1, components.sizeY - 1, 1, 2);
    }
}


template <class T>
T Matrix3<T>::trace() const
{
    if (sizeZ != 1)
        throw std::domain_error("Cannot compute the trace of the matrix : Matrix should be 2D (sizeZ = 1) but here sizeZ = " + std::to_string(this->sizeZ));
    if (sizeX != sizeY)
        throw std::domain_error("Cannot compute the trace of the matrix : Matrix should be squared (sizeX = sizeY) but here sizeX = " + std::to_string(this->sizeX) + " and sizeY = " + std::to_string(this->sizeY));
    T sum = T();
    for (size_t i = 0; i < sizeX; i++)
        sum += this->unsafe(i, i);
    return sum;
}

template <class T>
Matrix3<T>& Matrix3<T>::init(const std::vector<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ)
{
    this->data = data;
    this->sizeX = sizeX;
    this->sizeY = sizeY;
    this->sizeZ = sizeZ;
    _updateStrides();
    // std::cout << "Copied " << this->data.size() << std::endl;
    return *this;
}

template <class T>
std::ostream& operator<<(std::ostream& io, const Matrix3<T>& v) {
    io << v.toString();
    return io;
}

template <class T>
std::ostream& operator<<(std::ostream& io, std::shared_ptr<Matrix3<T>> v) {
    io << v->toString();
    return io;
}

template <class T>
T Matrix3<T>::min() const
{
    T min = std::numeric_limits<T>::max();
    for(const T& val : this->data)
        min = std::min(min, val);
    return min;
}
template <class T>
T Matrix3<T>::max() const
{
    T max = std::numeric_limits<T>::min();
    for(const T& val : this->data)
        max = std::max(max, val);
    return max;
}

template <class T>
Matrix3<T>& Matrix3<T>::min(const T &minVal)
{
    this->iterateParallel([&](size_t i) {
        this->data[i] = std::min(this->data[i], minVal);
    });
    return *this;
}

template <class T>
Matrix3<T>& Matrix3<T>::max(const T &maxVal)
{
    this->iterateParallel([&](size_t i) {
        this->data[i] = std::max(this->data[i], maxVal);
    });
    return *this;
}

template <class T>
Matrix3<T> Matrix3<T>::abs() const
{
    Matrix3<T> m = *this;
    for (T& val : m)
        val = std::abs(val);
    return m;
}

template <class T>
T Matrix3<T>::sum() const
{
    T sum = T();
    for (const auto& val : this->data)
        sum += val;
    return sum;
}

template <class T>
Matrix3<T> Matrix3<T>::rounded(int precision) const
{
    for(T& val : this->data)
        val = (val * std::pow(10, precision)) / (int)std::pow(10, precision);
    return *this;
}

template <class T>
Matrix3<T>& Matrix3<T>::normalize() {
    if (this->data.empty()) return *this;
    T min = data[0], max = data[0];
    for (const T& val : data)
    {
        if (min > val) min = val;
        if (max < val) max = val;
    }
    T diff = max - min;
    if (min != max) {
        for (T& val : data)
        {
            val = (val - min)/ diff;
        }
    } else {
        for (T& val : data)
        {
            val = T();
        }
    }
    return *this;
}
template <class T>
Matrix3<T> Matrix3<T>::normalized() const {
    Matrix3 mat = *this;
    return mat.normalize();
}

template <class T>
Matrix3<T>& Matrix3<T>::normalizeUsing(NORMALIZE_METHOD normalizeMethod) {
    if (this->data.empty()) return *this;
    if (normalizeMethod == NORMALIZE_MINMAX) {
        this->normalize();
    } else if (normalizeMethod == NORMALIZE_Z_SCORE) {
        auto [mu, sigma] = stats::getMuSigma(this->data);
        if (sigma == T()) return *this;
        *this = (*this - mu) / sigma;
    } else if (normalizeMethod == NORMALIZE_SOFTMAX) {
        Matrix3<T> exponentialSelf = *this;
        for (auto& val : exponentialSelf)
            val = std::exp(val);
        T exponentialSum = exponentialSelf.sum();
        if (exponentialSum == T()) return *this;
        *this = exponentialSelf / exponentialSum;
    }
    return *this;
}
template <class T>
Matrix3<T> Matrix3<T>::normalizedUsing(NORMALIZE_METHOD normalizeMethod) const {
    Matrix3 mat = *this;
    return mat.normalizeUsing(normalizeMethod);
}

template <class T>
Matrix3<T> Matrix3<T>::transposeXY()
{
    Matrix3<T> res(this->getDimensions().yxz());
    res.iterateParallel([&](int x, int y, int z) {
        res(x, y, z) = this->unsafe(y, x, z);
    });
    return res;
}

template <class T>
Matrix3<float> Matrix3<T>::gaussian(int sizeOnX, int sizeOnY, int sizeOnZ, float sigma, const Vector3& offset) {
    Matrix3<float> gaussian(sizeOnX, sizeOnY, sizeOnZ);
    Vector3 center = Vector3(sizeOnX/2.f, sizeOnY/2.f, sizeOnZ/2.f) + offset;
    center -= Vector3((sizeOnX > 1 ? .5 : 0), (sizeOnY > 1 ? .5 : 0), (sizeOnZ > 1 ? .5 : 0));
    float oneOverSqrt2Pi = 1.f/std::sqrt(2 * M_PI);
    float sqrSigma = sigma * sigma;
    gaussian.iterateParallel([&](const Vector3i& pos) {
        gaussian(pos) = std::exp(-((pos - center).norm2())/(2*sqrSigma)) * oneOverSqrt2Pi;
    });
    return gaussian;
}

template <class T>
Matrix3<float> Matrix3<T>::normalizedGaussian(int sizeOnX, int sizeOnY, int sizeOnZ, float sigma, const Vector3 &offset)
{
    Matrix3<float> gauss = Matrix3<float>::gaussian(sizeOnX, sizeOnY, sizeOnZ, sigma, offset);
    float sum = gauss.sum();
    if (sum != 0) gauss /= sum;
//    return gauss / (sum != 0 ? sum : 1.f);
    return gauss;
}

template <class T>
Matrix3<T> Matrix3<T>::LaplacianOfGaussian(int sizeOnX, int sizeOnY, int sizeOnZ, float sigma) const {
    Matrix3<float> laplacian = Matrix3<float>(3 + sizeOnX/2 +1, 3 + sizeOnY/2 +1, 3 + sizeOnZ/2 +1, 0.f);
    laplacian.at(sizeOnX/2 + 1, sizeOnY/2 + 1, sizeOnZ/2 + 1) = 1.f;
    laplacian = laplacian.laplacian();
    Matrix3<float> gaussian = Matrix3<T>::gaussian(sizeOnX, sizeOnY, sizeOnZ, sigma);
    return this->convolution(laplacian.convolution(gaussian));
}
template <class T>
Matrix3<T> Matrix3<T>::meanSmooth(int sizeOnX, int sizeOnY, int sizeOnZ, [[maybe_unused]] bool ignoreBorders) const {
    auto mirror_idx = [&](int i, size_t n) {
        if (n <= 1) return 0;
        // small-radius typical, so loop is tiny near borders
        while ((unsigned)i >= (unsigned)n) {
            if (i < 0) i = -i - 1;
            else       i = 2*n - i - 1;
        }
        return i;
    };

    const size_t X = this->sizeX;
    const size_t Y = this->sizeY;
    const size_t Z = this->sizeZ;

    const size_t rx = sizeOnX / 2;
    const size_t ry = sizeOnY / 2;
    const size_t rz = sizeOnZ / 2;

    // Scratch buffers
    Matrix3<T> buf1(X, Y, Z);
    Matrix3<T> buf2(X, Y, Z);

    const T* src = this->data.data();
    T* tmp  = buf1.data.data();
    T* tmp2 = buf2.data.data();

    const int XY = X * Y;

    // ---- Pass 1: X (contiguous) ----
#pragma omp parallel for collapse(2) schedule(static)
    for (size_t z = 0; z < Z; ++z) {
        for (size_t y = 0; y < Y; ++y) {
            const size_t base = z*XY + y*X;

            // initial window sum at x=0
            T sum = T();
            for (int dx = -int(rx); dx <= int(rx); ++dx) {
                int xi = mirror_idx(dx, X);
                sum += src[base + xi];
            }
            tmp[base + 0] = sum;

            // slide
            for (size_t x = 1; x < X; ++x) {
                int x_add = mirror_idx(x + rx, X);
                int x_sub = mirror_idx(x - rx - 1, X);
                sum += src[base + x_add] - src[base + x_sub];
                tmp[base + x] = sum;
            }
        }
    }

    // ---- Pass 2: Y (stride = X) ----
#pragma omp parallel for collapse(2) schedule(static)
    for (size_t z = 0; z < Z; ++z) {
        for (size_t x = 0; x < X; ++x) {
            const size_t baseZ = z*XY;

            // initial sum at y=0
            T sum = T();
            for (int dy = -int(ry); dy <= int(ry); ++dy) {
                int yi = mirror_idx(dy, Y);
                sum += tmp[baseZ + yi*X + x];
            }
            tmp2[baseZ + 0*X + x] = sum;

            for (size_t y = 1; y < Y; ++y) {
                int y_add = mirror_idx(y + ry, Y);
                int y_sub = mirror_idx(y - ry - 1, Y);
                sum += tmp[baseZ + y_add*X + x] - tmp[baseZ + y_sub*X + x];
                tmp2[baseZ + y*X + x] = sum;
            }
        }
    }

    // ---- Pass 3: Z (stride = X*Y) ----
    Matrix3<T> out(X, Y, Z);
    T* dst = out.data.data();


#pragma omp parallel for collapse(2) schedule(static)
    for (size_t y = 0; y < Y; ++y) {
        for (size_t x = 0; x < X; ++x) {

            // initial sum at z=0
            T sum = T();
            for (int dz = -int(rz); dz <= int(rz); ++dz) {
                int zi = mirror_idx(dz, Z);
                sum += tmp2[zi*XY + y*X + x];
            }
            dst[0*XY + y*X + x] = sum;

            for (size_t z = 1; z < Z; ++z) {
                int z_add = mirror_idx(z + rz, Z);
                int z_sub = mirror_idx(z - rz - 1, Z);
                sum += tmp2[z_add*XY + y*X + x] - tmp2[z_sub*XY + y*X + x];
                dst[z*XY + y*X + x] = sum;
            }
        }
    }
    const float inv = 1.0f / float(sizeOnX * sizeOnY * sizeOnZ);

    // scale in-place (or fuse into last pass if you prefer)
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < X*Y*Z; ++i) dst[i] *= inv;
    return out;
}

template <class T>
Matrix3<T> Matrix3<T>::gaussianSmooth(float sigma, bool ignoreZ, bool ignoreBorders, float limitFactor) const
{
    const size_t X = this->sizeX;
    const size_t Y = this->sizeY;
    const size_t Z = this->sizeZ;
    const size_t XY = X * Y;

    auto mirror_idx = [&](int i, int n) {
        if (n <= 1) return 0;
        while ((unsigned)i >= (unsigned)n) {
            if (i < 0) i = -i - 1;
            else       i = 2*n - i - 1;
        }
        return i;
    };

    if (sigma <= 0.0f) return *this;

    if (limitFactor < 0.0f) {
        limitFactor = std::max({X, Y, Z}) / sigma;
    }

    const size_t r = std::max(1, (int)std::ceil(limitFactor * sigma));
    const size_t K = 2*r + 1;

    // Build 1D kernel (centered at 0)
    std::vector<float> kernel(K);
    {
        const float inv2s2 = 1.0f / (2.0f * sigma * sigma);
        float sum = 0.0f;
        for (int k = -int(r); k <= int(r); ++k) {
            float w = std::exp(-(float)(k*k) * inv2s2);
            kernel[k + r] = w;
            sum += w;
        }
        const float invSum = 1.0f / sum;
        for (float& w : kernel) w *= invSum;
    }

    Matrix3<T> buf1(X, Y, Z);
    Matrix3<T> buf2(X, Y, Z);
    Matrix3<T> out (X, Y, Z);

    const T* src = this->data.data();
    T* tmp  = buf1.data.data();
    T* tmp2 = buf2.data.data();
    T* dst  = out.data.data();

// ---- Pass 1: X ----
#pragma omp parallel for collapse(2) schedule(static)
    for (size_t z = 0; z < Z; ++z) {
        for (size_t y = 0; y < Y; ++y) {
            const size_t base = z*XY + y*X;

            // Interior (no borders)
            for (size_t x = r; x < X - r; ++x) {
                T acc = T();
                const T* p = src + base + (x - r);
                for (size_t k = 0; k < K; ++k) {
                    acc += p[k] * kernel[k];
                }
                tmp[base + x] = (T)acc;
            }

            // Left border
            for (size_t x = 0; x < std::min(r, X); ++x) {
                T acc = T();
                for (int k = -int(r); k <= int(r); ++k) {
                    int xi = ignoreBorders ? (x + k) : mirror_idx(x + k, X);
                    if (ignoreBorders && ((unsigned)xi >= (unsigned)X)) continue;
                    acc += src[base + xi] * kernel[k + r];
                }
                tmp[base + x] = (T)acc;
            }

            // Right border
            for (size_t x = std::max(int(X - r), 0); x < X; ++x) {
                T acc = T();
                for (int k = -int(r); k <= int(r); ++k) {
                    int xi = ignoreBorders ? (x + k) : mirror_idx(x + k, X);
                    if (ignoreBorders && ((unsigned)xi >= (unsigned)X)) continue;
                    acc += src[base + xi] * kernel[k + r];
                }
                tmp[base + x] = (T)acc;
            }
        }
    }

// ---- Pass 2: Y ----
#pragma omp parallel for collapse(2) schedule(static)
    for (size_t z = 0; z < Z; ++z) {
        for (size_t x = 0; x < X; ++x) {
            const size_t baseZ = z*XY;

            // Interior
            for (size_t y = r; y < Y - r; ++y) {
                T acc = T();
                int idx = baseZ + (y - r)*X + x;
                for (size_t k = 0; k < K; ++k) {
                    acc += tmp[idx + k*X] * kernel[k];
                }
                tmp2[baseZ + y*X + x] = (T)acc;
            }

            // Borders
            for (size_t y = 0; y < std::min(r, Y); ++y) {
                T acc = T();
                for (int k = -int(r); k <= int(r); ++k) {
                    int yi = ignoreBorders ? (y + k) : mirror_idx(y + k, Y);
                    if (ignoreBorders && ((unsigned)yi >= (unsigned)Y)) continue;
                    acc += tmp[baseZ + yi*X + x] * kernel[k + r];
                }
                tmp2[baseZ + y*X + x] = (T)acc;
            }
            for (size_t y = std::max(int(Y - r), 0); y < Y; ++y) {
                T acc = T();
                for (int k = -int(r); k <= int(r); ++k) {
                    int yi = ignoreBorders ? (y + k) : mirror_idx(y + k, Y);
                    if (ignoreBorders && ((unsigned)yi >= (unsigned)Y)) continue;
                    acc += tmp[baseZ + yi*X + x] * kernel[k + r];
                }
                tmp2[baseZ + y*X + x] = (T)acc;
            }
        }
    }

    if (ignoreZ) {
        // If Z ignored, output is tmp2
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < X*Y*Z; ++i) dst[i] = tmp2[i];
        return out;
    }

// ---- Pass 3: Z ----
#pragma omp parallel for collapse(2) schedule(static)
    for (size_t y = 0; y < Y; ++y) {
        for (size_t x = 0; x < X; ++x) {

            // Interior
            for (size_t z = r; z < Z - r; ++z) {
                T acc = T();
                int idx = (z - r)*XY + y*X + x;
                for (size_t k = 0; k < K; ++k) {
                    acc += tmp2[idx + k*XY] * kernel[k];
                }
                dst[z*XY + y*X + x] = (T)acc;
            }

            // Borders
            for (size_t z = 0; z < std::min(r, Z); ++z) {
                T acc = T();
                for (int k = -int(r); k <= int(r); ++k) {
                    int zi = ignoreBorders ? (z + k) : mirror_idx(z + k, Z);
                    if (ignoreBorders && ((unsigned)zi >= (unsigned)Z)) continue;
                    acc += tmp2[zi*XY + y*X + x] * kernel[k + r];
                }
                dst[z*XY + y*X + x] = (T)acc;
            }
            for (size_t z = std::max(int(Z - r), 0); z < Z; ++z) {
                T acc = T();
                for (int k = -int(r); k <= int(r); ++k) {
                    int zi = ignoreBorders ? (z + k) : mirror_idx(z + k, Z);
                    if (ignoreBorders && ((unsigned)zi >= (unsigned)Z)) continue;
                    acc += tmp2[zi*XY + y*X + x] * kernel[k + r];
                }
                dst[z*XY + y*X + x] = (T)acc;
            }
        }
    }

    return out;
}

/*
template <class T>
Matrix3<T> Matrix3<T>::gaussianSmooth(float sigma, bool ignoreZ, bool ignoreBorders, float limitFactor) const
{
    ignoreBorders = false;
    int sizeX = this->sizeX;
    int sizeY = this->sizeY;
    int sizeZ = this->sizeZ;

    Matrix3<T> result(sizeX, sizeY, sizeZ);
    result.raiseErrorOnBadCoord = false;
    result.returned_value_on_outside = this->returned_value_on_outside;
    Matrix3<T> tempResult = *this; //(sizeX, sizeY, sizeZ);
    tempResult.raiseErrorOnBadCoord = false;
    tempResult.returned_value_on_outside = this->returned_value_on_outside;

    auto gaussian = [](float x, float sigma) {
        return exp(-(x * x) / (2 * sigma * sigma)) / (sqrt(2 * M_PI) * sigma);
    };

    float sumX = 0.0f;
    float sumY = 0.0f;
    float sumZ = 0.0f;

    if (limitFactor < 0)
        limitFactor = std::max({sizeX, sizeY, sizeZ}) / sigma;

    int startX = std::floor(sizeX / 2 - sigma * limitFactor);
    startX = (ignoreBorders ? startX : std::max(startX, 0));
    int endX = std::ceil(sizeX / 2 + sigma * limitFactor);
    endX = (ignoreBorders ? endX : std::min(endX, sizeX));

    std::vector<float> kernelX(sizeX);

    for (int i = 0; i < sizeX; ++i) {
        float x = i - sizeX / 2;
        if (std::abs(x) > limitFactor * sigma) continue;
        kernelX[i] = gaussian(x, sigma);
        sumX += kernelX[i];
    }

    // Convolve along x-axis
    iterateParallel([&](int x, int y, int z) {
        T sum = T();
        for (int dx = startX; dx < endX; ++dx) {
            int nx = x + dx - sizeX / 2;
            sum += tempResult.at(nx, y, z) * kernelX[dx];
        }
        result.at(x, y, z) = sum / sumX;
    });


    int startY = std::floor(sizeY / 2 - sigma * limitFactor);
    startY = (ignoreBorders ? startY : std::max(startY, 0));
    int endY = std::ceil(sizeY / 2 + sigma * limitFactor);
    endY = (ignoreBorders ? endY : std::min(endY, sizeY));

    std::vector<float> kernelY(sizeY);

    for (int i = 0; i < sizeY; ++i) {
        float x = i - sizeY / 2;
        if (std::abs(x) > limitFactor * sigma) continue;
        kernelY[i] = gaussian(x, sigma);
        sumY += kernelY[i];
    }

    // Convolve along y-axis
    result.iterateParallel([&](int x, int y, int z) {
        T sum = T();
        for (int dy = startY; dy < endY; ++dy) {
            int ny = y + dy - sizeY / 2;
            // if (ny >= 0 && ny < sizeY) {
                sum += result.at(x, ny, z) * kernelY[dy];
            // }
        }
        tempResult.at(x, y, z) = sum / sumY;
    });
    result = tempResult;

    if (!ignoreZ) {

        int startZ = std::floor(sizeZ / 2 - sigma * limitFactor);
        startZ = (ignoreBorders ? startZ : std::max(startZ, 0));
        int endZ = std::ceil(sizeZ / 2 + sigma * limitFactor);
        endZ = (ignoreBorders ? endZ : std::min(endZ, sizeZ));

        std::vector<float> kernelZ(sizeZ);

        for (int i = 0; i < sizeZ; ++i) {
            float x = i - sizeZ / 2;
            if (std::abs(x) > limitFactor * sigma) continue;
            kernelZ[i] = gaussian(x, sigma);
            sumZ += kernelZ[i];
        }

        // Convolve along z-axis
        tempResult.iterateParallel([&](int x, int y, int z) {
            T sum = T();
            for (int dz = startZ; dz < endZ; ++dz) {
                int nz = z + dz - sizeZ / 2;
                // if (nz >= 0 && nz < sizeZ) {
                    sum += tempResult.at(x, y, nz) * kernelZ[dz];
                // }
            }
            result.at(x, y, z) = sum / sumZ;
        });
    }

    return result;
}
*/

template <class T>
Matrix3<T> Matrix3<T>::medianBlur(int sizeOnX, int sizeOnY, int sizeOnZ, bool ignoreBorders) const
{
    Matrix3<T> result(this->getDimensions());

    int halfX = sizeOnX / 2;
    int halfY = sizeOnY / 2;
    int halfZ = sizeOnZ / 2;

    this->iterateParallel([&](const Vector3& p) {
        std::vector<T> values;
        values.reserve(sizeOnX * sizeOnY * sizeOnZ);
        for (int dx = -halfX; dx <= halfX; dx++) {
            for (int dy = -halfY; dy <= halfY; dy++) {
                for (int dz = -halfZ; dz <= halfZ; dz++) {
                    Vector3 pos(p.x() + dx, p.y() + dy, p.z() + dz);
                    if (!ignoreBorders || checkCoord(pos))
                        values.push_back(this->at(pos));
                }
            }
        }
        std::sort(values.begin(), values.end());
        if (values.size() % 2 == 1) {
            result(p) = values[values.size() / 2];
        } else {
            result(p) = (values[values.size() / 2 - 1] + values[values.size() / 2]) * .5f;
        }
    });
    return result;
}

template <class T>
Matrix3<Vector3> Matrix3<T>::gradient() const {
    Matrix3<Vector3> returningGrid(this->sizeX, this->sizeY, this->sizeZ);

    auto res = returningGrid.data.data();
    const auto src = this->data.data();

    const auto repeat = [&](int x, int y, int z) -> size_t {
        return this->checkCoord(x, y, z) ? getIndex(x, y, z) : this->getRepeatIndex(x, y, z);
    };
    if (this->sizeZ > 1) {
        iterateParallel([&](int x, int y, int z) {
            res[getIndex(x, y, z)] = Vector3((src[repeat(x + 1, y, z)] - src[repeat(x - 1, y, z)]) * .5f,
                                                (src[repeat(x, y + 1, z)] - src[repeat(x, y - 1, z)]) * .5f,
                                                (src[repeat(x, y, z + 1)] - src[repeat(x, y, z - 1)]) * .5f);
        });
    } else {
        iterateParallel([&](int x, int y, int z) {
            res[getIndex(x, y, z)] = Vector3((src[repeat(x + 1, y, z)] - src[repeat(x - 1, y, z)]) * .5f,
                                                (src[repeat(x, y + 1, z)] - src[repeat(x, y - 1, z)]) * .5f,
                                                0);
        });
    }
    return returningGrid;
}
template <class T>
inline Matrix3<Vector3> Matrix3<T>::grad() const {
    return this->gradient();
}

template <class T>
Matrix3<T> Matrix3<T>::laplacian() const
{
    Matrix3 returningGrid = *this;
    this->raiseErrorOnBadCoord = false;
    this->defaultValueOnBadCoord = T();
    iterateParallel([&](int x, int y, int z) {
        T val = T();
        val += this->at(x    , y    , z + 1);
        val += this->at(x    , y    , z - 1);
        val += this->at(x    , y + 1, z    );
        val += this->at(x    , y - 1, z    );
        val += this->at(x + 1, y    , z    );
        val += this->at(x - 1, y    , z    );
        val -= this->at(x    , y    , z    ) * 6;
        returningGrid.at(x, y, z) = val;
    });
    /*for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                T val = T();
                val += this->at(x    , y    , z + 1);
                val += this->at(x    , y    , z - 1);
                val += this->at(x    , y + 1, z    );
                val += this->at(x    , y - 1, z    );
                val += this->at(x + 1, y    , z    );
                val += this->at(x - 1, y    , z    );
                val -= this->at(x    , y    , z    ) * 6;
                returningGrid.at(x, y, z) = val;
            }
        }
    }*/
    this->raiseErrorOnBadCoord = true;
    return returningGrid;
}

template <class T>
Matrix3<int> Matrix3<T>::binarize(T limitValue, bool greaterValuesAreSetToOne, bool useAlsoTheEqualSign) const
{
    Matrix3<int> bin(this->sizeX, this->sizeY, this->sizeZ);
    iterateParallel([&](size_t i){
        if (greaterValuesAreSetToOne) {
            if ((useAlsoTheEqualSign && this->data[i] >= limitValue) || (!useAlsoTheEqualSign && this->data[i] > limitValue)) {
                bin[i] = 1;
            } else {
                bin[i] = 0;
            }
        } else {
            if ((useAlsoTheEqualSign && this->data[i] <= limitValue) || (!useAlsoTheEqualSign && this->data[i] < limitValue)) {
                bin[i] = 1;
            } else {
                bin[i] = 0;
            }
        }
    });
    /*for (size_t i = 0; i < this->size(); i++) {
        if (greaterValuesAreSetToOne) {
            if ((useAlsoTheEqualSign && this->data[i] >= limitValue) || (!useAlsoTheEqualSign && this->data[i] > limitValue)) {
                bin[i] = 1;
            } else {
                bin[i] = 0;
            }
        } else {
            if ((useAlsoTheEqualSign && this->data[i] <= limitValue) || (!useAlsoTheEqualSign && this->data[i] < limitValue)) {
                bin[i] = 1;
            } else {
                bin[i] = 0;
            }
        }
    }*/
    return bin;
}

template <class T>
Matrix3<int> Matrix3<T>::binarizeBetween(T minValue, T maxValue, bool insideValuesAreSetToOne, bool useAlsoTheEqualSign) const
{
    Matrix3<int> bin(this->sizeX, this->sizeY, this->sizeZ);
    iterateParallel([&](size_t i){
        if (insideValuesAreSetToOne) {
            if ((useAlsoTheEqualSign && this->data[i] >= minValue && this->data[i] <= maxValue) || (!useAlsoTheEqualSign && this->data[i] > minValue && this->data[i] < maxValue)) {
                bin[i] = 1;
            } else {
                bin[i] = 0;
            }
        } else {
            if ((useAlsoTheEqualSign && this->data[i] <= minValue && this->data[i] <= maxValue) || (!useAlsoTheEqualSign && this->data[i] < minValue && this->data[i] < maxValue)) {
                bin[i] = 1;
            } else {
                bin[i] = 0;
            }
        }
    });

    /*for (size_t i = 0; i < this->size(); i++) {
        if (insideValuesAreSetToOne) {
            if ((useAlsoTheEqualSign && this->data[i] >= minValue && this->data[i] <= maxValue) || (!useAlsoTheEqualSign && this->data[i] > minValue && this->data[i] < maxValue)) {
                bin[i] = 1;
            } else {
                bin[i] = 0;
            }
        } else {
            if ((useAlsoTheEqualSign && this->data[i] <= minValue && this->data[i] <= maxValue) || (!useAlsoTheEqualSign && this->data[i] < minValue && this->data[i] < maxValue)) {
                bin[i] = 1;
            } else {
                bin[i] = 0;
            }
        }
    }*/
    return bin;
}

template <class T>
Matrix3<int> Matrix3<T>::isosurface(T isovalue, bool ignoreZtopBorder, bool ignoreBorders) const
{
    Matrix3<int> surface = Matrix3<int>(this->getDimensions());
    Matrix3<T> copy = *this;
    copy.raiseErrorOnBadCoord = false;
    copy.defaultValueOnBadCoord = T();
    bool useZ = (this->sizeZ > 1);

    iterateParallel([&](int x, int y, int z) {
        if (this->unsafe(x, y, z) <= 0) return;
        if (ignoreBorders && (x == 0 || x == int(sizeX) - 1 || y == 0 || y == int(sizeY) - 1 || z == 0 || (ignoreZtopBorder && z == int(sizeZ) - 1))) return;
        bool isSurface = false;
        for (int dx = -1; dx <= 1 && !isSurface; dx++) {
            for (int dy = -1; dy <= 1 && !isSurface; dy++) {
                for (int dz = (useZ ? -1 : 0); dz <= (useZ ? 1 : 0); dz++) {
                    if ((this->at(x + dx, y + dy, z + dz) - isovalue) <= 0)
                        isSurface = true;
                }
            }
        }
        surface.at(x, y, z) = (isSurface ? 1 : 0);
    });
    return surface;
}

template <class T>
Matrix3<T> Matrix3<T>::slice(int index, int axis) const
{
    Matrix3<T> result;
    if (axis == 0) { // YZ
        result = Matrix3<T>(1, sizeY, sizeZ);
        for (size_t y = 0; y < sizeY; y++) {
            for (size_t z = 0; z < sizeZ; z++) {
                result(0, y, z) = this->unsafe(index, y, z);
            }
        }
    } else if (axis == 1) { // XZ
        result = Matrix3<T>(sizeX, 1, sizeZ);
        for (size_t x = 0; x < sizeX; x++) {
            for (size_t z = 0; z < sizeZ; z++) {
                result(x, 0, z) = this->unsafe(x, index, z);
            }
        }
    } else if (axis == 2) { // XY
        result = Matrix3<T>(sizeX, sizeY, 1);
        for (size_t x = 0; x < sizeX; x++) {
            for (size_t y = 0; y < sizeY; y++) {
                result(x, y, 0) = this->unsafe(x, y, index);
            }
        }
    }
    return result;
}

template <class T>
Matrix3<T> Matrix3<T>::sliceXY(int index) const
{
    return slice(index, 2);
}

template <class T>
Matrix3<T> Matrix3<T>::sliceYZ(int index) const
{
    return slice(index, 0);
}

template <class T>
Matrix3<T> Matrix3<T>::sliceXZ(int index) const
{
    return slice(index, 1);
}

/*template <class T>
Matrix3<T> Matrix3<T>::random(const Vector3& dimensions)
{
    return Matrix3<T>::random(dimensions.x(), dimensions.y(), dimensions.z());
}*/
template <class T>
Matrix3<T> Matrix3<T>::random(const Vector3i& dimensions)
{
    return Matrix3<T>::random(dimensions.x(), dimensions.y(), dimensions.z());
}

template <class T>
Matrix3<T> Matrix3<T>::random(size_t sizeX, size_t sizeY, size_t sizeZ)
{
    Matrix3<T> res(sizeX, sizeY, sizeZ);
    res.iterateParallel([&](size_t i) {
        res[i] = random_gen::generate();
    });
    /*for (auto& v : res)
        v = random_gen::generate();*/
    return res;
}

template <class T>
Matrix3<T>& Matrix3<T>::insertRow(size_t indexToInsert, int affectedDimension, T newData)
{
    auto it = this->data.begin();
    int jumps = 0;
    switch(affectedDimension) {
    case 0:
        jumps = sizeX;
        if (indexToInsert < 0) indexToInsert = sizeX; // If default value, it's last X-index
        it += indexToInsert;
        for (; it <= this->data.end() - (indexToInsert == 0 ? 1 : 0); it += jumps) {
            it = this->data.insert(it, newData) + 1; // Set "it" to the value next to the inserted value
        }
        this->sizeX++;
        break;
    case 1:
        jumps = sizeX * sizeY;
        if (indexToInsert < 0) indexToInsert = sizeY; // If default value, it's last Y-index
        it += (indexToInsert * sizeX);
        for (; it <= this->data.end() - (indexToInsert == 0 ? 1 : 0); it += jumps) {
            for (size_t i = 0; i < sizeX; i++) {
                it = this->data.insert(it, newData) + 1;
            }
        }
        this->sizeY++;
        break;
    case 2:
        if (indexToInsert < 0) indexToInsert = sizeZ; // If default value, it's last Z-index
        it += (sizeX * sizeY) * indexToInsert;
        for (size_t i = 0; i < sizeX * sizeY; i++)
            it = this->data.insert(it, newData) + 1;
        this->sizeZ++;
        break;
    default:
        throw std::out_of_range("insertRow can only be processed on dim 0, 1 or 2 (resp. X, Y, Z)");
    }
    return *this;
}

template <class T>
Matrix3<T> Matrix3<T>::identity(size_t sizeX, size_t sizeY, size_t sizeZ)
{
    static_assert(std::is_arithmetic<T>::value, "");
    Matrix3<T> mat(sizeX, sizeY, sizeZ);
    if (sizeZ == 1) {
        if (sizeX != sizeY)
            throw std::invalid_argument("Identity matrix must be square (dim X == dim Y)");
        for (size_t i = 0; i < sizeX; i++)
            mat.at(i, i, 0) = 1;
    } else {
        if (sizeX != sizeY || sizeX != sizeZ)
            throw std::invalid_argument("All dimensions must be the same length to do a 3-d identity matrix");
        for (size_t i = 0; i < sizeX; i++)
            mat.at(i, i, i) = 1;
    }
    return mat;
}

/*template <class T>
Matrix3<T> Matrix3<T>::perlin(const Vector3 &dimensions, const Vector3& scale, int seed)
{
    Matrix3<T> result(dimensions);

    result.iterateParallel([&](float x, float y, float z) {
        result(x, y, z) = random_gen::generate_perlin(x * scale.x() + seed, y * scale.y() + seed, z * scale.z() + seed);
    });
    return result;
}*/
template <class T>
Matrix3<T> Matrix3<T>::perlin(const Vector3i& dimensions, const Vector3& scale, int seed)
{
    Matrix3<T> result(dimensions);

    result.iterateParallel([&](float x, float y, float z) {
        result(x, y, z) = random_gen::generate_perlin(x * scale.x() + seed, y * scale.y() + seed, z * scale.z() + seed);
    });
    return result;
}

template <class T>
Matrix3<T> operator+(Matrix3<T> a, const Matrix3<T> &b) {
    a += b;
    return a;
}
template <class T>
Matrix3<T> Matrix3<T>::operator-() const {
    return *this * -1.f;
}

template <class T> template <class U>
Matrix3<T>& Matrix3<T>::operator+=(const Matrix3<U>& o) {
    if (this->sizeX != o.sizeX || this->sizeY != o.sizeY || this->sizeZ != o.sizeZ)
        throw std::domain_error("Matrices must have same sizes to be added (M1 = " + this->toString() + " and M2 = " + o.toString());
    iterateParallel([&](size_t i) {
        data[i] += o.data[i];
    });
//    for (size_t i = 0; i < data.size(); i++) {
//        data[i] += o.data[i];
//    }
    return *this;
}
template <class T, typename U>
Matrix3<T> operator-(Matrix3<T> a, const Matrix3<U> &b) {
    a -= b;
    return a;
}
template <class T> template <class U>
Matrix3<T>& Matrix3<T>::operator-=(const Matrix3<U>& o)  {
    if (this->sizeX != o.sizeX || this->sizeY != o.sizeY || this->sizeZ != o.sizeZ)
        throw std::domain_error("Matrices must have same sizes to be substracted (M1 = " + this->toString() + " and M2 = " + o.toString());
    iterateParallel([&](size_t i) {
        data[i] -= o.data[i];
    });
    /*for (size_t i = 0; i < data.size(); i++) {
        data[i] -= o.data[i];
    }*/
    return *this;
}
template <class T, typename U>
Matrix3<T> operator*(Matrix3<T> a, const Matrix3<U>& o) {
    a *= o;
    return a;
}
template <class T> template <class U>
Matrix3<T>& Matrix3<T>::operator*=(const Matrix3<U> &o) {
    if (this->sizeX != o.sizeX || this->sizeY != o.sizeY || this->sizeZ != o.sizeZ)
        throw std::domain_error("Matrices must have same sizes to be multiplied (M1 = " + this->toString() + " and M2 = " + o.toString());
    iterateParallel([&](size_t i) {
        data[i] *= o.data[i];
    });
    /*
    for (size_t i = 0; i < data.size(); i++) {
        data[i] *= o.data[i];
    }*/
    return *this;
}
template <class T, typename U>
Matrix3<T> operator/(Matrix3<T> a, const Matrix3<U>& b) {
    a /= b;
    return a;
}
template <class T> template <class U>
Matrix3<T>& Matrix3<T>::operator/=(const Matrix3<U>& o) {
    if (this->sizeX != o.sizeX || this->sizeY != o.sizeY || this->sizeZ != o.sizeZ)
        throw std::domain_error("Matrices must have same sizes to be divided (M1 = " + this->toString() + " and M2 = " + o.toString());
    iterateParallel([&](size_t i) {
        data[i] /= o.data[i];
    });
    /*
    for (size_t i = 0; i < data.size(); i++) {
        data[i] /= o.data[i];
    }*/
    return *this;
}
template <class T, typename U>
Matrix3<T> operator*(Matrix3<T> a, U o) {
    a *= o;
    return a;
}
template <class T> template <class U>
Matrix3<T>& Matrix3<T>::operator*=(U o) {
    iterateParallel([&](size_t i) {
        data[i] *= o;
    });
    /*
    for (size_t i = 0; i < data.size(); i++) {
        data[i] *= o;
    }*/
    return *this;
}

template <class T, typename U>
Matrix3<T> operator/(Matrix3<T> a, U o) {
    a /= o;
    return a;
}
template <class T> template <class U>
Matrix3<T>& Matrix3<T>::operator/=(U o) {
    iterateParallel([&](size_t i) {
        data[i] /= o;
    });
    /*
    for (size_t i = 0; i < data.size(); i++) {
        data[i] /= o;
    }*/
    return *this;
}
template <class T, typename U>
Matrix3<T> operator+(Matrix3<T> a, U o) {
    a += o;
    return a;
}
template <class T> template <class U>
Matrix3<T>& Matrix3<T>::operator+=(U o) {
    iterateParallel([&](size_t i) {
        data[i] += o;
    });
    /*
    for (size_t i = 0; i < data.size(); i++) {
        data[i] += o;
    }*/
    return *this;
}
template <class T, typename U>
Matrix3<T> operator-(Matrix3<T> a, U o) {
    a -= o;
    return a;
}
template <class T> template <class U>
Matrix3<T>& Matrix3<T>::operator-=(U o) {
    iterateParallel([&](size_t i) {
        data[i] -= o;
    });
    /*
    for (size_t i = 0; i < data.size(); i++) {
        data[i] -= o;
    }*/
    return *this;
}
template <class T> template <class U>
bool Matrix3<T>::operator==(const Matrix3<U> &o) const {
    if (this->sizeX != o.sizeX || this->sizeY != o.sizeY || this->sizeZ != o.sizeZ)
        return false;
    for (size_t i = 0; i < this->data.size(); i++)
        if (this->data[i] != o.data[i])
            return false;
    return true;
}

template <class T>
int Matrix3<T>::getNumberNeighbors(size_t x, size_t y, size_t z, bool using4connect) const
{
    int neighbors = 0;
    if (using4connect) {
        if (x > 0) neighbors++;
        if (x < sizeX-1) neighbors++;
        if (y > 0) neighbors++;
        if (y < sizeY - 1) neighbors++;
        if (z > 0) neighbors++;
        if (z < sizeZ - 1) neighbors++;
    }
    else {
        for (int dx = -1; dx <= 1; dx++)
            for (int dy = -1; dy <= 1; dy++)
                for (int dz = -1; dz <= 1; dz++)
                    if (dx != 0 || dy != 0 || dz != 0)
                        neighbors += (checkCoord(x + dx, y + dy, z + dz) ? 1 : 0);
    }
    return neighbors;
}
/*template <class T>
int Matrix3<T>::getNumberNeighbors(const Vector3& pos, bool using4connect) const
{
    return getNumberNeighbors(pos.x(), pos.y(), pos.z(), using4connect);
}*/
template <class T>
int Matrix3<T>::getNumberNeighbors(const Vector3i& pos, bool using4connect) const
{
    return getNumberNeighbors(pos.x(), pos.y(), pos.z(), using4connect);
}

template <class T>
Matrix3<T> Matrix3<T>::resize(float factor, RESIZE_MODE mode) const
{
    Vector3 newSize = this->getDimensions() * factor;
    if (newSize.x() < 1) newSize.x() = 1;
    if (newSize.y() < 1) newSize.y() = 1;
    if (newSize.z() < 1) newSize.z() = 1;
    return this->resize(newSize, mode);
}
/*template <class T>
Matrix3<T> Matrix3<T>::resize(const Vector3& newSize, RESIZE_MODE mode) const
{
    return this->resize(newSize.x(), newSize.y(), newSize.z(), mode);
}*/
template <class T>
Matrix3<T> Matrix3<T>::resize(const Vector3i& newSize, RESIZE_MODE mode) const
{
    return this->resize(newSize.x(), newSize.y(), newSize.z(), mode);
}

template <class T>
Matrix3<T> Matrix3<T>::resize(size_t newX, size_t newY, size_t newZ, RESIZE_MODE mode) const
{
    Matrix3<T> newMat(newX, newY, newZ);
    newMat.raiseErrorOnBadCoord = false;
    float rx = (this->sizeX - 1) / std::max(1.f, (float)(newX - 1)), ry = (this->sizeY - 1) / std::max(1.f, (float)(newY - 1)), rz = (this->sizeZ - 1) / std::max(1.f, (float)(newZ - 1));

    auto dst = newMat.data.data();
    //const auto src = data.data();

    if (mode == LINEAR) {
        newMat.iterateParallel([&](size_t x, size_t y, size_t z) {
            // float d_x = (float(x) * rx) - int(x * rx);
            // float d_y = (float(y) * ry) - int(y * ry);
            // float d_z = (float(z) * rz) - int(z * rz);
            dst[newMat.getIndex(x, y, z)] = interpolate(x * rx, y * ry, z * rz);
            /*size_t x_original = size_t(x * rx);
            size_t x_plus_1 = (x_original >= this->sizeX - 1 ? x_original : x_original + 1);
            float d_x = (x * rx) - x_original;
            size_t y_original = size_t(y * ry);
            size_t y_plus_1 = (y_original >= this->sizeY - 1 ? y_original : y_original + 1);
            float d_y = (y * ry) - y_original;
            size_t z_original = size_t(z * rz);
            size_t z_plus_1 = (z_original >= this->sizeZ - 1 ? z_original : z_original + 1);
            float d_z = (z * rz) - z_original;

            T f000 = this->unsafe(x_original    , y_original    , z_original    );
            T f100 = this->unsafe(x_plus_1, y_original    , z_original    );
            T f010 = this->unsafe(x_original    , y_plus_1, z_original    );
            T f110 = this->unsafe(x_plus_1, y_plus_1, z_original    );
            T f001 = this->unsafe(x_original    , y_original    , z_plus_1);
            T f101 = this->unsafe(x_plus_1, y_original    , z_plus_1);
            T f011 = this->unsafe(x_original    , y_plus_1, z_plus_1);
            T f111 = this->unsafe(x_plus_1, y_plus_1, z_plus_1);
            // Interpolation
            T res = ((
                                      f000 * (1-d_x) + f100 * d_x) * (1-d_y) + (
                                      f010 * (1-d_x) + f110 * d_x) * d_y) * (1 - d_z) +
                                ((
                                     f001 * (1-d_x) + f101 * d_x) * (1-d_y) + (
                                     f011 * (1-d_x) + f111 * d_x) * d_y) * d_z;

            newMat.at(x, y, z) = res;*/
        });
    } else if (mode == NEAREST) {
        newMat = this->resizeNearest(newX, newY, newZ);

    } else if (mode == MAX_VAL || mode == MIN_VAL) {
//        for (auto& val : newMat)
//            val = (mode == MAX_VAL ? std::numeric_limits<T>::min() : std::numeric_limits<T>::max());
        Matrix3<short int> modifiedMatrix(newX, newY, newZ, 0);
        iterateParallel([&](int x, int y, int z) {
            size_t startX = x / rx;
            size_t endX = (x + 1) / rx;
            size_t startY = y / ry;
            size_t endY = (y + 1) / ry;
            size_t startZ = z / rz;
            size_t endZ = (z + 1) / rz;

            // Not sure that this is the most efficient strategy, but meh...
            for (size_t dx = startX; dx <= endX; dx++) {
                for (size_t dy = startY; dy <= endY; dy++) {
                    for (size_t dz = startZ; dz <= endZ; dz++) {
                        if (modifiedMatrix.checkCoord(dx, dy, dz)) {
                            // If this cell hasn't been modified yet, we cannot apply the min/max operator
                            if (modifiedMatrix.at(dx, dy, dz) == 0) {
                                newMat.at(dx, dy, dz) = this->unsafe(x, y, z);
                                modifiedMatrix.at(dx, dy, dz) = 1;
                            } else {
                                newMat.at(dx, dy, dz) = (mode == MAX_VAL ? std::max(newMat.at(dx, dy, dz), this->unsafe(x, y, z)) : std::min(newMat.at(dx, dy, dz), this->unsafe(x, y, z)));
                            }
                        }
                    }
                }
            }
        });
    } else if (mode == FILL_WITH_DEFAULT) {
        if (rz == 0) rz = 1;
        Vector3 ratio(rx, ry, rz);
        iterateParallel([&](const Vector3& pos) {
            newMat(pos / ratio) = this->unsafe(pos);
        });
    }
    newMat.raiseErrorOnBadCoord = this->raiseErrorOnBadCoord;
    return newMat;
}

template <class T>
Matrix3<T> Matrix3<T>::resizeNearest(float factor) const
{
    Vector3 newSize = this->getDimensions() * factor;
    if (newSize.x() < 1) newSize.x() = 1;
    if (newSize.y() < 1) newSize.y() = 1;
    if (newSize.z() < 1) newSize.z() = 1;
    return this->resizeNearest(newSize);
}

template <class T>
Matrix3<T> Matrix3<T>::resizeNearest(size_t newX, size_t newY, size_t newZ) const
{
    Matrix3<T> newMat(newX, newY, newZ);
    newMat.raiseErrorOnBadCoord = false;
    float rx = (this->sizeX - 1) / std::max(1.f, (float)(newX - 1)), ry = (this->sizeY - 1) / std::max(1.f, (float)(newY - 1)), rz = (this->sizeZ - 1) / std::max(1.f, (float)(newZ - 1));

    // Apply interpolations
    Vector3 ratio(rx, ry, rz);
    newMat.iterateParallel([&](const Vector3i& pos) {
        newMat(pos) = this->unsafe((pos * ratio).roundedDown());
    });
    return newMat;
}

template <class T>
Matrix3<T> Matrix3<T>::resizeNearest(const Vector3i& newSize) const
{
    return this->resizeNearest(newSize.x(), newSize.y(), newSize.z());
}

template <class T>
Matrix3<T> Matrix3<T>::subset(const Vector3i& start, const Vector3i& end) const
{
    int endZ = end.z();
    if (start.z() == 0 && end.z() == 0)
        endZ = -1; // Give it the default value so it will be managed by the main function
    return this->subset(start.x(), end.x(), start.y(), end.y(), start.z(), endZ);
}

template <class T>
Matrix3<T> Matrix3<T>::subset(int startX, int endX, int startY, int endY, int startZ, int endZ) const
{
    if (endZ == -1) endZ = this->sizeZ;
    if (startX == 0 && endX == sizeX && startY == 0 && endY == sizeY && startZ == 0 &&  endZ == sizeZ) return *this;

    Matrix3<T> croppedMatrix(std::max(endX - startX, 0), std::max(endY - startY, 0), std::max(endZ - startZ, 0));
    croppedMatrix.iterateParallel([&](int x, int y, int z) {
        int oldX = x + startX;
        int oldY = y + startY;
        int oldZ = z + startZ;
        if (0 > oldX || oldX >= this->sizeX || 0 > oldY || oldY >= this->sizeY || 0 > oldZ || oldZ >= this->sizeZ) return;
        croppedMatrix(x, y, z) = this->unsafe(oldX, oldY, oldZ);
    });
    return croppedMatrix;
}

template <class T>
inline Matrix3<T>& Matrix3<T>::paste(const Matrix3<T> &matrixToPaste, const Vector3i& upperLeftFrontCorner)
{
    return this->paste(matrixToPaste, upperLeftFrontCorner.x(), upperLeftFrontCorner.y(), upperLeftFrontCorner.z());
}
template <class T>
Matrix3<T>& Matrix3<T>::paste(const Matrix3<T>& matrixToPaste, int left, int up, int front)
{
    auto dst = data.data();
    const auto src = matrixToPaste.data.data();

    iterateParallel([&](int x, int y, int z) {
       int oldX = x - left;
       int oldY = y - up;
       int oldZ = z - front;
       if (!checkCoord(x, y, z) || !matrixToPaste.checkCoord(oldX, oldY, oldZ)) return;
       dst[getIndex(x, y, z)] = src[matrixToPaste.getIndex(oldX, oldY, oldZ)];
    });
    return *this;
}

template <class T>
Matrix3<T>& Matrix3<T>::add(const Matrix3<T>& matrixToAdd, const Vector3& upperLeftFrontCorner, [[maybe_unused]] bool useInterpolation)
{
    if (useInterpolation) {
        matrixToAdd.iterate([&](const Vector3i& pos) {
            this->addValueAt(matrixToAdd(pos), upperLeftFrontCorner + pos);
        });
        return *this;
    } else {
        return this->add(matrixToAdd, upperLeftFrontCorner.x(), upperLeftFrontCorner.y(), upperLeftFrontCorner.z(), useInterpolation);
    }
}
template <class T>
Matrix3<T>& Matrix3<T>::add(const Matrix3<T> &matrixToAdd, int left, int up, int front, [[maybe_unused]] bool useInterpolation)
{
    auto dst = data.data();
    const auto src = matrixToAdd.data.data();

    matrixToAdd.iterateParallel([&](int x, int y, int z) {
       int oldX = x + left;
       int oldY = y + up;
       int oldZ = z + front;
       if (!checkCoord(oldX, oldY, oldZ) || !matrixToAdd.checkCoord(x, y, z)) return;
       dst[getIndex(oldX, oldY, oldZ)] += src[matrixToAdd.getIndex(x, y, z)];
    });
    return *this;
}

template <class T>
Matrix3<T> Matrix3<T>::concat(const Matrix3<T>& matrixToConcat)
{
    Matrix3<T> newMatrix(this->getDimensions() + matrixToConcat.getDimensions() * Vector3(1, 0, 0));
    newMatrix.paste(*this, Vector3());
    newMatrix.paste(matrixToConcat, (newMatrix.getDimensions() - matrixToConcat.getDimensions()) * Vector3(1, 0, 0));
    return newMatrix;
}

template <class T>
inline Matrix3<T>& Matrix3<T>::max(const Matrix3<T>& otherMatrix, const Vector3i& upperLeftFrontCorner)
{
    return this->max(otherMatrix, upperLeftFrontCorner.x(), upperLeftFrontCorner.y(), upperLeftFrontCorner.z());
}
template <class T>
Matrix3<T>& Matrix3<T>::max(const Matrix3<T>& otherMatrix, int left, int up, int front)
{
    auto dst = data.data();
    const auto src = otherMatrix.data.data();

    otherMatrix.iterateParallel([&](int x, int y, int z) {
       int oldX = x + left;
       int oldY = y + up;
       int oldZ = z + front;
       if (!checkCoord(oldX, oldY, oldZ) || !otherMatrix.checkCoord(x, y, z)) return;
       size_t idx = getIndex(oldX, oldY, oldZ);
       dst[idx] = std::max(dst[idx], src[otherMatrix.getIndex(x, y, z)]);
    });
    return *this;
}

template <class T>
inline Matrix3<T>& Matrix3<T>::min(const Matrix3<T>& otherMatrix, const Vector3i& upperLeftFrontCorner)
{
    return this->min(otherMatrix, upperLeftFrontCorner.x(), upperLeftFrontCorner.y(), upperLeftFrontCorner.z());
}
template <class T>
Matrix3<T>& Matrix3<T>::min(const Matrix3<T> &otherMatrix, int left, int up, int front)
{
    auto dst = data.data();
    const auto src = otherMatrix.data.data();

    otherMatrix.iterateParallel([&](int x, int y, int z) {
       int oldX = x + left;
       int oldY = y + up;
       int oldZ = z + front;
       if (!checkCoord(oldX, oldY, oldZ) || !otherMatrix.checkCoord(x, y, z)) return;
       size_t idx = getIndex(oldX, oldY, oldZ);
       dst[idx] = std::min(dst[idx], src[otherMatrix.getIndex(x, y, z)]);
    });
    return *this;
}

template <class T>
Matrix3<T> Matrix3<T>::max(const Matrix3<T>& m1, const Matrix3<T>& m2)
{
    if (m1.getDimensions() != m2.getDimensions())
        throw std::domain_error("Matrices must have same sizes to be maxed (M1 = " + m1.toString() + " and M2 = " + m2.toString());
    Matrix3<T> res(m1.getDimensions());
    auto R = res.data.data();
    const auto M1 = m1.data.data();
    const auto M2 = m2.data.data();
    res.iterateParallel([&](size_t i) {
        R[i] = std::max(M1[i], M2[i]);
    });
    return res;
}

template <class T>
Matrix3<T> Matrix3<T>::min(const Matrix3<T>& m1, const Matrix3<T>& m2)
{
    if (m1.getDimensions() != m2.getDimensions())
        throw std::domain_error("Matrices must have same sizes to be mined (M1 = " + m1.toString() + " and M2 = " + m2.toString());
    Matrix3<T> res(m1.getDimensions());
    auto R = res.data.data();
    const auto M1 = m1.data.data();
    const auto M2 = m2.data.data();
    res.iterateParallel([&](size_t i) {
        R[i] = std::min(M1[i], M2[i]);
    });
    return res;
}

template <class T>
Matrix3<float> Matrix3<T>::toDistanceMap(bool ignoreZlayer, bool considerBorders)
{
    Matrix3<float> distances(this->sizeX, this->sizeY, this->sizeZ, std::numeric_limits<float>::max() - 10000);
    distances.raiseErrorOnBadCoord = false;
    distances.defaultValueOnBadCoord = (considerBorders ? 0.f : std::numeric_limits<float>::max() - 10000);

    // Using the Chamfer distance -> direct neighbor               => distance = 3
    //                               diagonal on 2 axis neighbor   => distance = 4
    //                               diagonal on all axis neighbor => distance = 5
//    float predefinedDistances[4] = {0, 3, 4, 5};
    // First pass
    distances.iterate([&](const Vector3& pos) {
        float currentVal = distances.at(pos);
        if (!this->unsafe(pos)) {
            distances.at(pos) = 0;
            return;
        }
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                for (int dz = -1; dz <= 1; dz++) {
                    if (ignoreZlayer && dz != 0) continue;
                    // Weighted distance transform
//                            currentVal = std::min(currentVal, distances.at(dx, dy, dz) + predefinedDistances[std::abs(dx) + std::abs(dy) + std::abs(dz)]);
                    currentVal = std::min(currentVal, distances.at(pos.x()+dx, pos.y()+dy, pos.z()+dz) + (float)std::sqrt(dx*dx + dy*dy + dz*dz));
                }
            }
        }
        distances.at(pos) = currentVal;
    });
    // Second pass
    distances.iterateReverse([&](const Vector3& pos) {
        if (!this->unsafe(pos)) {
            distances.at(pos) = 0;
            return;
        }
        float currentVal = distances.at(pos);
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                for (int dz = -1; dz <= 1; dz++) {
                    if (ignoreZlayer && dz != 0) continue;
                    currentVal = std::min(currentVal, distances.at(pos.x()+dx, pos.y()+dy, pos.z()+dz) + (float)std::sqrt(dx*dx + dy*dy + dz*dz));
                }
            }
        }
        distances.at(pos) = currentVal;
    });
    distances.raiseErrorOnBadCoord = true;
    distances.defaultValueOnBadCoord = 0.f;
    return distances; //.normalize();
}

template <class T>
Matrix3<std::complex<float>> Matrix3<T>::FFT() const {
    int dimX = (isPowerOf2(sizeX) ? sizeX : findNextPowerOfTwo(sizeX));
    int dimY = (isPowerOf2(sizeY) ? sizeY : findNextPowerOfTwo(sizeY));
    int dimZ = (isPowerOf2(sizeZ) ? sizeZ : findNextPowerOfTwo(sizeZ));
    Matrix3<std::complex<float>> result0(dimX, dimY, dimZ); // Result for X-axis FFT
    Matrix3<std::complex<float>> result1(dimX, dimY, dimZ); // Result for Y-axis FFT

    auto R0 = result0.data.data();
    auto R1 = result1.data.data();

    const size_t strideY = result0._strideY;
    const size_t strideZ = result0._strideZ;
    auto indx = [=](size_t x, size_t y, size_t z) { return z * strideZ + y * strideY + x; };

    // Perform FFT along X-axis for each YZ plane
    #pragma omp parallel for collapse(2)
    for (size_t y = 0; y < dimY; ++y) {
        for (size_t z = 0; z < dimZ; ++z) {
            std::vector<std::complex<float>> col(dimX);
            for (size_t x = 0; x < dimX; ++x) {
                col[x] = std::complex<float>((checkCoord(x, y, z) ? float(this->unsafe(x, y, z)) : 0.f), 0);
            }
            std::vector<std::complex<float>> fft_result_x = fft(col); // Perform 1D FFT along X-axis
            for (size_t x = 0; x < dimX; ++x) {
                R0[indx(x, y, z)] = fft_result_x[x];
            }
        }
    }

    // Perform FFT along Y-axis for each XZ plane
    #pragma omp parallel for collapse(2)
    for (size_t x = 0; x < dimX; ++x) {
        for (size_t z = 0; z < dimZ; ++z) {
            std::vector<std::complex<float>> col(dimY);
            for (size_t y = 0; y < dimY; ++y) {
                col[y] = R0[indx(x, y, z)]; // Use the result from X-axis FFT
            }
            std::vector<std::complex<float>> fft_result_y = fft(col); // Perform 1D FFT along Y-axis
            for (size_t y = 0; y < dimY; ++y) {
                R1[indx(x, y, z)] = fft_result_y[y];
            }
        }
    }

    // Perform FFT along Z-axis for each XY plane
    #pragma omp parallel for collapse(2)
    for (size_t x = 0; x < dimX; ++x) {
        for (size_t y = 0; y < dimY; ++y) {
            std::vector<std::complex<float>> row(dimZ);
            for (size_t z = 0; z < dimZ; ++z) {
                row[z] = R1[indx(x, y, z)]; // Use the result from Y-axis FFT
            }
            std::vector<std::complex<float>> fft_result_z = fft(row); // Perform 1D FFT along Z-axis
            for (size_t z = 0; z < dimZ; ++z) {
                R0[indx(x, y, z)] = fft_result_z[z];
            }
        }
    }

    return result0; // Return the final result after X, Y, and Z axis FFTs
}
template <class T>
Matrix3<std::complex<float> > Matrix3<T>::iFFT(const Vector3i& finalDimensions) const
{
    if (!finalDimensions.isValid()) {
        return this->iFFT(sizeX, sizeY, sizeZ);
    } else {
        return this->iFFT(finalDimensions.x(), finalDimensions.y(), finalDimensions.z());
    }
}
template <class T>
Matrix3<std::complex<float> > Matrix3<T>::iFFT(size_t cropX, size_t cropY, size_t cropZ) const
{
    auto initial = Matrix3<std::complex<float>>(findNextPowerOfTwo(sizeX), findNextPowerOfTwo(sizeY), findNextPowerOfTwo(sizeZ));
    initial.paste(*this, 0, 0, 0);
    Matrix3<std::complex<float>> result(this->getDimensions()); // Result after inverse FFT

    auto R = result.data.data();

    // Perform inverse FFT along X-axis for each YZ plane
    #pragma omp parallel for collapse(2)
    for (size_t y = 0; y < initial.sizeY; ++y) {
        for (size_t z = 0; z < initial.sizeZ; ++z) {
            std::vector<std::complex<float>> col(initial.sizeX);
            for (size_t x = 0; x < initial.sizeX; ++x) {
                col[x] = initial.unsafe(x, y, z);
            }
            std::vector<std::complex<float>> ifft_result_x = inverseFFT(col); // Perform 1D IFFT along X-axis
            for (size_t x = 0; x < initial.sizeX; ++x) {
                R[result.getIndex(x, y, z)] = ifft_result_x[x];
            }
        }
    }

    // Perform inverse FFT along Y-axis for each XZ plane
    #pragma omp parallel for collapse(2)
    for (size_t x = 0; x < initial.sizeX; ++x) {
        for (size_t z = 0; z < initial.sizeZ; ++z) {
            std::vector<std::complex<float>> col(initial.sizeY);
            for (size_t y = 0; y < initial.sizeY; ++y) {
                col[y] = R[result.getIndex(x, y, z)]; // Use the result after X-axis inverse FFT
            }
            std::vector<std::complex<float>> ifft_result_y = inverseFFT(col); // Perform 1D IFFT along Y-axis
            for (size_t y = 0; y < initial.sizeY; ++y) {
                R[result.getIndex(x, y, z)] = ifft_result_y[y];
            }
        }
    }

    // Perform inverse FFT along Z-axis for each XY plane
    #pragma omp parallel for collapse(2)
    for (size_t x = 0; x < initial.sizeX; ++x) {
        for (size_t y = 0; y < initial.sizeY; ++y) {
            std::vector<std::complex<float>> row(initial.sizeZ);
            for (size_t z = 0; z < initial.sizeZ; ++z) {
                row[z] = R[result.getIndex(x, y, z)]; // Use the result after Y-axis inverse FFT
            }
            std::vector<std::complex<float>> ifft_result_z = inverseFFT(row); // Perform 1D IFFT along Z-axis
            for (size_t z = 0; z < initial.sizeZ; ++z) {
                R[result.getIndex(x, y, z)] = ifft_result_z[z];
            }
        }
    }

    return result.subset(0, cropX, 0, cropY, 0, cropZ); // Return the final result after X, Y, and Z axis inverse FFTs
}

template <class T>
Matrix3<T> Matrix3<T>::flip(bool onX, bool onY, bool onZ)
{
    Matrix3<T> result = *this;

    auto dst = result.data.data();
    const auto src = data.data();
    result.iterateParallel([&](int x, int y, int z) {
        int targetX = (onX ? this->sizeX - (x +1) : x);
        int targetY = (onY ? this->sizeY - (y +1) : y);
        int targetZ = (onZ ? this->sizeZ - (z +1) : z);
        dst[result.getIndex(x, y, z)] = src[getIndex(targetX, targetY, targetZ)];
    });
    return result;
}

template <class T> template <class U>
Matrix3<T> Matrix3<T>::convolution(const Matrix3<U>& convMatrix, CONVOLUTION_BORDERS borders) const
{
    //return convolutionIgnoredBorders(*this, convMatrix);
    Matrix3<T> result(this->sizeX, this->sizeY, this->sizeZ);
//    this->raiseErrorOnBadCoord = false;

    // Pre-calculate normalisation value
    U normalisationValue = convMatrix.sum();

    // Choose border handling method once before loop
    auto handleBorder = [&](Vector3& pos) {
        if (borders == CONVOLUTION_BORDERS::IGNORED && !result.checkCoord(pos))
            return false;
        if (borders == CONVOLUTION_BORDERS::MIRROR && !result.checkCoord(pos))
            pos = getMirrorPosition(pos);
        else if (borders == CONVOLUTION_BORDERS::REPEAT && !result.checkCoord(pos))
            pos = getRepeatPosition(pos);
        else if (borders == CONVOLUTION_BORDERS::WRAPPING && !result.checkCoord(pos))
            pos = this->getWrappedPosition(pos);
        return true;
    };

    auto getVal = [&](const Vector3& pos) {
        if (borders == CONVOLUTION_BORDERS::ZERO_PAD && !result.checkCoord(pos))
            return T();
        return this->unsafe(pos);
    };

    this->iterateParallel([&](int x, int y, int z) {
        U divisor = normalisationValue;
        T neighboringSum = T();
        convMatrix.iterate([&](int dx, int dy, int dz) {
            int dt_x = dx - (convMatrix.sizeX / 2);
            int dt_y = dy - (convMatrix.sizeY / 2);
            int dt_z = dz - (convMatrix.sizeZ / 2);
            Vector3 cellValuePosition(x + dt_x, y + dt_y, z + dt_z);

            if (handleBorder(cellValuePosition)) {
                neighboringSum += (handleBorder(cellValuePosition) ? getVal(cellValuePosition) * convMatrix.at(dx, dy, dz) : T());
            } else {
                divisor -= convMatrix(dx, dy, dz);
            }
        });
        result.at(x, y, z) = neighboringSum;
        if (divisor != U())
            result.at(x, y, z) /= divisor;
    });
    return result;
}

template <class T>
inline Vector3 Matrix3<T>::getMirrorPosition(const Vector3& pos)  const noexcept
{
    float x = pos.x();
    float y = pos.y();
    float z = pos.z();
    x = int(x < 0 ? std::abs(x) : (x >= sizeX ? (sizeX - 1) - (x - sizeX) : x));
    y = int(y < 0 ? std::abs(y) : (y >= sizeY ? (sizeY - 1) - (y - sizeY) : y));
    z = int(z < 0 ? std::abs(z) : (z >= sizeZ ? (sizeZ - 1) - (z - sizeZ) : z));
    return Vector3(x, y, z);
}

template <class T>
inline Vector3 Matrix3<T>::getWrappedPosition(const Vector3& pos) const noexcept
{
    Vector3 rounded = pos.roundedDown();
    Vector3 decimals = pos - rounded;
    Vector3  wrap = Vector3(int(rounded.x() + sizeX) % sizeX,
                           int(rounded.y() + sizeY) % sizeY,
                           int(rounded.z() + sizeZ) % sizeZ
                           ) + decimals;
    return wrap;
}

template <class T>
inline Vector3 Matrix3<T>::getRepeatPosition(const Vector3& pos) const noexcept
{
    Vector3 returned;
    returned.x() = std::min(std::max(0.f, pos.x()), (float)sizeX - 1);
    returned.y() = std::min(std::max(0.f, pos.y()), (float)sizeY - 1);
    returned.z() = std::min(std::max(0.f, pos.z()), (float)sizeZ - 1);
    return returned;
}

template <class T>
inline Vector3i Matrix3<T>::getMirrorPosition(const Vector3i& pos)  const noexcept
{
    int x = pos.x();
    int y = pos.y();
    int z = pos.z();
    x = (x < 0 ? std::abs(x) : (x >= int(sizeX) ? (int(sizeX) - 1) - (x - sizeX) : x));
    y = (y < 0 ? std::abs(y) : (y >= int(sizeY) ? (int(sizeY) - 1) - (y - sizeY) : y));
    z = (z < 0 ? std::abs(z) : (z >= int(sizeZ) ? (int(sizeZ) - 1) - (z - sizeZ) : z));
    return Vector3i(x, y, z);
}

template<class T>
inline size_t Matrix3<T>::getMirrorIndex(int x, int y, int z) const noexcept
{
    return getIndex(x < 0 ? std::abs(x) : (x >= int(sizeX) ? (int(sizeX) - 1) - (x - sizeX) : x),
                    y < 0 ? std::abs(y) : (y >= int(sizeY) ? (int(sizeY) - 1) - (y - sizeY) : y),
                    z < 0 ? std::abs(z) : (z >= int(sizeZ) ? (int(sizeZ) - 1) - (z - sizeZ) : z));
}

template <class T>
inline Vector3i Matrix3<T>::getWrappedPosition(const Vector3i& pos) const noexcept
{
    return Vector3i((pos.x() + sizeX) % sizeX, (pos.y() + sizeY) % sizeY, (pos.z() + sizeZ) % sizeZ);
}

template<class T>
inline size_t Matrix3<T>::getWrappedIndex(int x, int y, int z) const noexcept
{
    return getIndex((x + sizeX) % sizeX, (y + sizeY) % sizeY, (z + sizeZ) % sizeZ);
}

template <class T>
inline Vector3i Matrix3<T>::getRepeatPosition(const Vector3i& pos) const noexcept
{
    Vector3i returned;
    returned.x() = std::min(std::max(0, pos.x()), int(sizeX) - 1);
    returned.y() = std::min(std::max(0, pos.y()), int(sizeY) - 1);
    returned.z() = std::min(std::max(0, pos.z()), int(sizeZ) - 1);
    return returned;
}

template<class T>
inline size_t Matrix3<T>::getRepeatIndex(int x, int y, int z) const noexcept
{
    return getIndex(std::min(std::max(0, x), int(sizeX) - 1), std::min(std::max(0, y), int(sizeY) - 1), std::min(std::max(0, z), int(sizeZ) - 1));
}

template <class T>
Matrix3<T> Matrix3<T>::warpWith(const Matrix3<Vector3>& warper, int precision) const
{
    // Warp definition : f(warp(p)) = f(p + warp vec)
    // But f(p) != f(p - warp vec), in the definition. I think it should create
    // better results because we can fetch outside values (with mirror for ex)
    // But that's not the current definition.
    Matrix3<T> result(getDimensions());
    // auto self = *this;
    // self.raiseErrorOnBadCoord = false;
    // self.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::DEFAULT_VALUE;
    // result.raiseErrorOnBadCoord = false;

    // auto valOrDefault = (this->returned_value_on_outside == RETURN_VALUE_ON_OUTSIDE::DEFAULT_VALUE ?
    //                      [&](const Vector3& q) -> T { return this->checkCoord(q) ? this->unsafe(q) : this->defaultValueOnBadCoord; } :
    //                      [&](const Vector3& q) -> T { return this->unsafe(this->getMirrorPosition(q)); });


    auto valOrDefault = (this->returned_value_on_outside == RETURN_VALUE_ON_OUTSIDE::DEFAULT_VALUE ?
                         std::function<T(const Vector3&)>([&](const Vector3& q) -> T { return this->checkCoord(q) ? this->unsafe(q) : this->defaultValueOnBadCoord; }) :
                         std::function<T(const Vector3&)>([&](const Vector3& q) -> T { return this->unsafe(this->getMirrorPosition(q)); }));

    if (precision > 1) {
        Matrix3<Vector3> displacements(warper.getDimensions());
        displacements.iterateParallel([&](const Vector3i& p) {
            auto& pos = displacements(p);
            pos = p;
            for (int iter = 0; iter < precision; iter++) {
                pos -= warper.interpolate(pos) / float(precision);
            }
        });
        iterateParallel([&](const Vector3& pos) {
            const auto& warp = displacements.at(pos);
            result.at(pos) = valOrDefault(warp);
        });
    } else {
        iterateParallel([&](const Vector3& pos) {
            const Vector3& warp = warper.at(pos);
            result.addValueAt(valOrDefault(pos), pos + warp);
        });
    }
    return result;
}

template <class T>
Matrix3<T> Matrix3<T>::warpWith(const CatmullRomSpline& original, const CatmullRomSpline& warperCurve) const
{
    // For now, start from a straight line on the X-axis
//    BSpline original = BSpline({this->getDimensions() * Vector3(0, .5, .5) + Vector3(1, 0, 0), this->getDimensions() * Vector3(1, .5, .5) - Vector3(1, 0, 0)});
    float pathsResolution = 1000.f;
    std::vector<Vector3> originalCurvePoints = original.getPath(pathsResolution);
    std::vector<Vector3> warperCurvePoints = warperCurve.getPath(pathsResolution);

    Matrix3<Vector3> warper(this->getDimensions());
    warper.raiseErrorOnBadCoord = false;
    Matrix3<float> modifications(this->getDimensions(), 0.f);
    modifications.raiseErrorOnBadCoord = false;

    // Vectors along the curve
    for (size_t i = 0; i < originalCurvePoints.size(); i++) {
        Vector3 pos = originalCurvePoints[i];
        Vector3 dir = warperCurvePoints[i] - pos;
        float curveWarpLength = dir.norm();
        dir.normalize();

        // In direction of the curve warping
        Vector3 endingPropagationPoint = Collision::intersectionRayAABBox(pos + Vector3(.5, .5, .5), dir, Vector3(), getDimensions());
        float distanceToBorder = (dir.norm2() > 0 ? (endingPropagationPoint - pos).norm() : 1.f);

        for (int j = 0; j < distanceToBorder; j++) {
            float warpLength = (1 - interpolation::linear(j / distanceToBorder)) * curveWarpLength;
//            warper.at(pos + dir * (float)j) += dir * warpLength;
//            modifications.at(pos + dir * (float)j) += 1.f;
            warper.addValueAt(dir * warpLength, pos + dir * (float)j);
            modifications.addValueAt(1.f, pos + dir * (float)j);
        }

        // In opposite direction of the curve warping
        endingPropagationPoint = Collision::intersectionRayAABBox(pos + Vector3(.5, .5, .5), dir * -1.f, Vector3(), getDimensions());
        distanceToBorder = (dir.norm2() > 0 ? (endingPropagationPoint - pos).norm() : 1.f);

        for (int j = 1; j < distanceToBorder; j++) {
            float warpLength = (1 - interpolation::linear(j / distanceToBorder)) * curveWarpLength;
//            warper.at(pos - dir * (float)j) += dir * warpLength;
//            modifications.at(pos - dir * (float)j) += 1.f;
            warper.addValueAt(dir * warpLength, pos - dir * (float)j);
            modifications.addValueAt(1.f, pos - dir * (float)j);
        }
    }

    for (size_t i = 0; i < warper.size(); i++) {
        warper[i] /= (modifications[i] > 0 ? modifications[i] : 1.f);
        modifications[i] = (modifications[i] > 0 ? 1.f : 0.f);
    }

    // Fill the empty gaps left
    while (modifications.min() <= 0) {
        for (size_t i = 0; i < warper.size(); i++) {
            if (modifications[i] == 0) {
                Vector3 pos = warper.getCoordAsVector3(i);
                Vector3 replaceValue;
                float divisor = 0.f;
                for (int x = -1; x <= 1; x++) {
                    for (int y = -1; y <= 1; y++) {
                        for (int z = -1; z <= 1; z++) {
                            Vector3 offset(x, y, z);
                            if (modifications.at(pos + offset) > 0) {
                                replaceValue += warper.at(pos + offset);
                                divisor++;
                            }
                        }
                    }
                }
                if (divisor > 0.f) {
                    warper.at(pos) = replaceValue / divisor;
                    modifications.at(pos) = 1.f;
                }

            }
        }
    }
//    Matrix3<float> tempReturn(warper.getDimensions());
//    for (size_t i = 0; i < tempReturn.size(); i++) {
//        tempReturn[i] = warper[i].x;
//    }
//    return tempReturn;
    return this->warpWith(warper);

    //    Matrix3<float> unit_ctrl(this->getDimensions(), 1.f);
}

template <class T>
Matrix3<T> Matrix3<T>::warpWithoutInterpolation(const Matrix3<Vector3>& warper) const
{
    Matrix3<T> result = *this; //(getDimensions());
    result.raiseErrorOnBadCoord = false;
    // auto self = *this;
    // self.raiseErrorOnBadCoord = false;
    // self.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;

    auto repeat = [&](const Vector3& q) -> Vector3 {
        return this->checkCoord(q) ? q : this->getRepeatPosition(q);
    };
    iterateParallel([&](const Vector3& pos) {
        result.at(pos + warper(pos)) = this->unsafe(repeat(pos));
    });
    return result;
}

template <class T>
Matrix3<T> Matrix3<T>::warpWithoutInterpolation(const CatmullRomSpline &original, const CatmullRomSpline &warperCurve) const
{
//    bool previousRaise = this->raiseErrorOnBadCoord;
//    T previousDefault  = this->returned_value_on_outside;
    // auto self = *this;
    // self.raiseErrorOnBadCoord = false;
    // self.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;

    auto repeat = [&](const Vector3& q) -> Vector3 {
        return this->checkCoord(q) ? q : this->getRepeatPosition(q);
    };

    Matrix3<Vector3> indices(this->getDimensions());
    for (size_t i = 0; i < this->size(); i++)
        indices[i] = indices.getCoordAsVector3(i);

    indices = indices.warpWith(original, warperCurve);

    Matrix3<T> values(this->getDimensions());
    for (size_t i = 0; i < this->size(); i++) {
        values[i] = this->unsafe(repeat(indices[i]));
    }

//    this->raiseErrorOnBadCoord = previousRaise;
//    this->returned_value_on_outside = previousDefault;

    return values;
}

template <class T>
Matrix3<float> Matrix3<T>::fbmNoise1D(FastNoiseLite noise, int sizeX, int sizeY, int sizeZ)
{
    Matrix3<float> values(sizeX, sizeY, sizeZ);
    values.iterateParallel([&](float x, float y, float z) {
        values(x, y, z) = noise.GetNoise(x, y, z);
    });
    return values;
}

template <class T>
Matrix3<Vector3> Matrix3<T>::fbmNoise2D(FastNoiseLite noise, int sizeX, int sizeY, int sizeZ)
{
    Matrix3<Vector3> values = Matrix3<Vector3>::fbmNoise3D(noise, sizeX, sizeY, sizeZ);
    values.iterateParallel([&] (size_t i) {
        values[i] = values[i].xy();
    });
    return values;
}

template <class T>
Matrix3<Vector3> Matrix3<T>::fbmNoise3D(FastNoiseLite noise, int sizeX, int sizeY, int sizeZ)
{
    Vector3 offsetDim1 = Vector3(   0,    0,    0);
    Vector3 offsetDim2 = Vector3(  42,  103, 2048);
    Vector3 offsetDim3 = Vector3(  15,  128, 1000);
    Matrix3<Vector3> values(sizeX, sizeY, sizeZ);
    values.iterateParallel([&](int x, int y, int z) {
        values(x, y, z) = Vector3(noise.GetNoise(x + offsetDim1.x(), y + offsetDim1.y(), z + offsetDim1.z()),
                               noise.GetNoise(x + offsetDim2.x(), y + offsetDim2.y(), z + offsetDim2.z()),
                               noise.GetNoise(x + offsetDim3.x(), y + offsetDim3.y(), z + offsetDim3.z()));
    });
    return values;
}

template <class T>
int Matrix3<T>::width() const {
    return this->sizeX;
}
template <class T>
int Matrix3<T>::depth() const {
    return this->sizeY;
}
template <class T>
int Matrix3<T>::height() const {
    return this->sizeZ;
}

template <class T>
int Matrix3<T>::rows() const
{
    return this->sizeY;
}
template <class T>
int Matrix3<T>::cols() const
{
    return this->sizeX;
}

template <class T>
Vector3i Matrix3<T>::getDimensions() const
{
    return Vector3i(this->sizeX, this->sizeY, this->sizeZ);
}


typedef Matrix3<float> GridF;
typedef Matrix3<int> GridI;
typedef Matrix3<Vector3> GridV3;


std::string stringifyGridF(const GridF& data, bool binaryMode = true);
std::string stringifyGridV3(const GridV3& data, bool binaryMode = true);

GridF loadGridF(const std::string& str, bool binaryMode = true);
GridV3 loadGridV3(const std::string& str, bool binaryMode = true);

#endif // MATRIX3_H
