#include "OMP_algo.h"

#include "Utils/Utils.h"
#include "DataStructure/Matrix3.h"
#include "DataStructure/Image.h"

#include "GUIElements/ImageViewer.h"

inline double dot(const GridF& a, const GridF& b) {
    return (a * b).sum();
}
inline double norm2(const GridF& a) {
    return std::sqrt(dot(a, a));
}

OMP::OMP()
{}

GridF OMP::getLargeReconstruction(const GridF &input, const Vector3i &finalDimensions)
{
    int scaleX = input.sizeX / tileSize;
    int scaleY = input.sizeY / tileSize;

    // Running OMP process on each tile
    std::vector<std::pair<Vector3, OMPResult>> results(scaleX * scaleY  * overlap  * overlap);
    #pragma omp parallel for collapse(2)
    for (int _i = 0; _i < scaleX  * overlap; _i++) {
        for (int _j = 0; _j < scaleY  * overlap; _j++) {
            float i = _i / float(overlap);
            float j = _j / float(overlap);
            auto targetPatch = input.subset(Vector3(tileSize * i, tileSize * j), Vector3(tileSize * (i + 1), tileSize * (j + 1)));
            OMPResult res = orthogonalMatchingPursuit(targetPatch, this->smallDictionary, nbPrimitives, 1e-5f);
            results[_i * scaleY * overlap + _j] = {Vector3(i, j), res};
        }
    }

    // Merge as a single image
    auto img = GridF(tileSize * scaleX * upscaleFactor, tileSize * scaleY * upscaleFactor, 1);
    GridF weights = GridF(img.getDimensions());
    GridF mask = GridF(tileSize * upscaleFactor, tileSize * upscaleFactor, 1, 1.f);
    Vector3 maskCenter = mask.getDimensions().xy() / 2.f;
    float maxRadius = maskCenter.norm();
    mask.iterateParallel([&](const Vector3& p) { mask[p] = std::clamp(1.f - interpolation::sigmoid((p - maskCenter).norm() / maxRadius), 0.01f, 1.f); });
    for (int iOverlap = 0; iOverlap < overlap; iOverlap++) {
        for (int jOverlap = 0; jOverlap < overlap; jOverlap++) {
            #pragma omp parallel for collapse(2)
            for (int x = 0; x < scaleX; x++) {
                for (int y = 0; y < scaleY; y++) {
                    int idx = (x * overlap + iOverlap) * scaleY * overlap + (y * overlap + jOverlap);
                    auto& [pos, res] = results[idx];
                    // std::cout << idx << ": " << pos << " -> (" << x << " * " << overlap << " + " << iOverlap << ") * " << scaleY << " * " << overlap << " + (" << y << " * " << overlap << " + " << jOverlap << ")" << std::endl;
                    GridF reconstruction = reconstructImage(res, this->bigDictionary);
                    img.add(reconstruction * mask, pos * tileSize * upscaleFactor);
                    weights.add(mask, pos * tileSize * upscaleFactor);
                }
            }
        }
    }
    img /= weights;
    if (finalDimensions != Vector3i::invalid && finalDimensions != img.getDimensions())
        return img.resize(finalDimensions);
    return img;
}


std::vector<double> OMP::solveLinearSystem(std::vector<std::vector<double>> M, std::vector<double> b)
{
    const int n = static_cast<int>(b.size());

    for (int k = 0; k < n; ++k)
    {
        // Pivot
        int pivot = k;
        double maxAbs = std::fabs(M[k][k]);
        for (int i = k + 1; i < n; ++i)
        {
            double v = std::fabs(M[i][k]);
            if (v > maxAbs)
            {
                maxAbs = v;
                pivot = i;
            }
        }

        if (maxAbs < 1e-12f)
            throw std::runtime_error("Singular system in OMP least-squares");

        if (pivot != k)
        {
            std::swap(M[k], M[pivot]);
            std::swap(b[k], b[pivot]);
        }

        // Elimination
        for (int i = k + 1; i < n; ++i)
        {
            double factor = M[i][k] / M[k][k];
            M[i][k] = 0.0f;
            for (int j = k + 1; j < n; ++j)
                M[i][j] -= factor * M[k][j];
            b[i] -= factor * b[k];
        }
    }

    // Back substitution
    std::vector<double> x(n, 0.0f);
    for (int i = n - 1; i >= 0; --i)
    {
        double sum = b[i];
        for (int j = i + 1; j < n; ++j)
            sum -= M[i][j] * x[j];
        x[i] = sum / M[i][i];
    }

    return x;
}

GridF OMP::combineAtoms(const std::vector<GridF>& dictionary, const std::vector<int>& support, const std::vector<double>& coeffs)
{
    if (support.empty())
        throw std::runtime_error("Support is empty");

    GridF out = dictionary[support[0]] * 0.0f;
    for (size_t i = 0; i < support.size(); ++i)
        out += dictionary[support[i]] * coeffs[i];
    return out;
}

std::vector<double> OMP::solveLeastSquaresActiveSet(const GridF& y, const std::vector<GridF>& dictionary, const std::vector<int>& support)
{
    const int m = static_cast<int>(support.size());

    std::vector<std::vector<double>> G(m, std::vector<double>(m, 0.0f));
    std::vector<double> rhs(m, 0.0f);
    #pragma omp parallel for
    for (int i = 0; i < m; ++i)
    {
        const GridF& ai = dictionary[support[i]];
        rhs[i] = dot(ai, y);
        #pragma omp parallel for
        for (int j = 0; j < m; ++j)
        {
            const GridF& aj = dictionary[support[j]];
            G[i][j] = dot(ai, aj);
        }
    }

    return solveLinearSystem(G, rhs);
}

OMPResult OMP::orthogonalMatchingPursuit(const GridF& y, const std::vector<GridF>& dictionary, int maxAtoms, double residualTolerance, bool debug)
{
    if (dictionary.empty())
        throw std::runtime_error("Dictionary is empty");
    if (maxAtoms <= 0)
        throw std::runtime_error("maxAtoms must be > 0");

    OMPResult result;
    result.residual = y;

    std::vector<bool> selected(dictionary.size(), false);

    for (int iter = 0; iter < maxAtoms; ++iter)
    {
        // 1) Select atom with maximum absolute correlation
        int bestIdx = -1;
        double bestCorr = -1.0f;

        for (size_t j = 0; j < dictionary.size(); ++j)
        {
            if (selected[j])
                continue;

            auto normalizedResidual = result.residual;
            double corr = std::fabs(dot(dictionary[j], normalizedResidual));
            if (corr > bestCorr)
            {
                bestCorr = corr;
                bestIdx = static_cast<int>(j);
            }
        }

        if (bestIdx < 0)
            break;

        result.support.push_back(bestIdx);
        selected[bestIdx] = true;

        // 2) Solve least-squares on current support
        result.coeffs = solveLeastSquaresActiveSet(y, dictionary, result.support);

        // 3) Update approximation and residual
        result.approximation = combineAtoms(dictionary, result.support, result.coeffs);
        result.residual = y - result.approximation;

        // 4) Stop if residual small enough
        if (dot(result.residual, result.residual) <= residualTolerance * residualTolerance)
            break;
    }

    return result;
}

/*
void OMP::createDictionary(int width, int height, DictionaryType type)
{
    std::vector<GridF> dict;

    switch (type)
    {
    case DictionaryType::Dirac:
        dict = createDictionary_Dirac(width, height);
        break;

    case DictionaryType::DCT:
        dict = createDictionary_DCT(width, height);
        break;

    case DictionaryType::Gaussian:
        dict = createDictionary_Gaussian(width, height, {1.0f, 2.0f, 4.0f});
        break;
    }

    normalizeDictionary(dict);
    return dict;
}

std::vector<GridF> OMP::createDictionary_Dirac(int width, int height)
{
    std::vector<GridF> dict(width * height);

    #pragma omp parallel for collapse(2)
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            GridF atom(width, height);
            atom.reset(0.0f);

            atom(x, y) = 1.0f;

            dict[y * width + x] = atom;
        }
    }
    return dict;
}

std::vector<GridF> OMP::createDictionary_DCT(int width, int height)
{
    std::vector<GridF> dict(width * height);

    #pragma omp parallel for collapse(2)
    for (int u = 0; u < width; ++u)
    {
        for (int v = 0; v < height; ++v)
        {
            GridF atom(width, height);

            double alpha_u = (u == 0) ? std::sqrt(1.0f / width) : std::sqrt(2.0f / width);
            double alpha_v = (v == 0) ? std::sqrt(1.0f / height) : std::sqrt(2.0f / height);

            atom.iterateParallel([&](size_t x, size_t y, size_t) {
                double value = alpha_u * alpha_v * std::cos((PI * (2 * x + 1) * u) / (2 * width)) * std::cos((PI * (2 * y + 1) * v) / (2 * height));

                atom.unsafe(x, y) = value;
            });

            dict[u * height + v] = atom;
        }
    }
    return dict;
}
std::vector<GridF> OMP::createDictionary_Gaussian(int width, int height, const std::vector<double>& sigmas)
{
    std::vector<GridF> dict(sigmas.size() * height * width);

    #pragma omp parallel for collapse(3)
    for (size_t iSigma = 0; iSigma < sigmas.size(); iSigma++)
    {
        for (int cy = 0; cy < height; ++cy)
        {
            for (int cx = 0; cx < width; ++cx)
            {
                auto& sigma = sigmas[iSigma];
                GridF atom(width, height);
                atom.iterateParallel([&](size_t x, size_t y, size_t) {
                    double dx = x - cx;
                    double dy = y - cy;

                    double val = std::exp(-(dx * dx + dy * dy) / (2.0f * sigma * sigma));
                    atom.unsafe(x, y) = val;
                });

                // normalize
                double n = norm2(atom);
                if (n > 1e-6f)
                    atom *= (1.0f / n);

                dict[iSigma * height * width + cx * height + cy] = atom;
            }
        }
    }
    return dict;
}
*/
std::vector<GridF> OMP::createDictionaryFromImage(int width, int height, const GridF &image)
{
    size_t nX = image.sizeX / width;
    size_t nY = image.sizeY / height;
    std::vector<GridF> patches(nX * nY, GridF(width, height));
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < nX; i++) {
        for (size_t j = 0; j < nY; j++) {
            GridF patch = image.subset(i * width, (i + 1) * width, j * height, (j + 1) * height, 0, 1);
            patches[i * nY + j].data = patch.data;
        }
    }
    normalizeDictionary(patches);
    return patches;
}

void OMP::createDictionaryFromImages(const std::vector<std::string>& paths)
{
    for (auto path : paths) {
        GridF bigImage = Image::readFromFile(path).getBwImage().resize(Vector3i(400, 400, 1));
        auto atoms = createDictionaryFromImage(tileSize * upscaleFactor, tileSize * upscaleFactor, bigImage);
        bigDictionary.insert(bigDictionary.end(), atoms.begin(), atoms.end());
    }
    std::shuffle(bigDictionary.begin(), bigDictionary.end(), random_gen::random_generator);
    bigDictionary.resize(std::min((int)bigDictionary.size(), maxAtoms));
    smallDictionary = bigDictionary;
    for (auto& data : smallDictionary) {
        data = data.resize(Vector3i(tileSize, tileSize, 1));
    }
}

void OMP::normalizeDictionary(std::vector<GridF>& dict)
{
    for (auto& atom : dict)
    {
        double n = norm2(atom);
        if (n > 1e-12f)
            atom *= (1.0f / n);
    }
}

GridF OMP::reconstructImage(const OMPResult& result, const std::vector<GridF>& dictionary)
{
    if (result.support.empty())
        throw std::runtime_error("OMPResult support is empty");

    if (result.support.size() != result.coeffs.size())
        throw std::runtime_error("support and coeffs size mismatch");

    GridF image = dictionary[result.support[0]] * 0.0f;

    for (size_t k = 0; k < result.support.size(); ++k)
    {
        int atomIndex = result.support[k];
        double coeff = result.coeffs[k];

        image += dictionary[atomIndex] * coeff;
    }

    return image;
}


/*
void testingSmoothmaxInCPP() {
    #define FLOAT_TYPE float
    int nbTrials = 10000;

    auto LSE = [](FLOAT_TYPE a, FLOAT_TYPE b, FLOAT_TYPE k) {
        return (1.f/k) * log(exp(k * a) + exp(k * b));
    };
    auto AbsMax = [](FLOAT_TYPE a, FLOAT_TYPE b, FLOAT_TYPE k) {
        return ((a + b)/2.f) + log(1.f + exp(k * abs(a - b)))/(2.f * k);
    };
    auto smaxLong = [](FLOAT_TYPE a, FLOAT_TYPE b, FLOAT_TYPE k) {
        return (a != b ? a + ((b - a)/2.f) * (1.f/(1 + exp(-k * (b - a))) + 1.f/(1 - exp(-k * (b - a)))) : a + 1.f/(2.f * k));
    };
    auto smaxMini = [](FLOAT_TYPE a, FLOAT_TYPE b, FLOAT_TYPE k) {
        return (a != b ? a + (b - a)/(1.f-exp(-2.f * k * (b - a))) : a + 1.f/(2.f * k));
    };

    FLOAT_TYPE a = 256.f;
    FLOAT_TYPE b = 0.f;
    FLOAT_TYPE k = 2.f;

    FLOAT_TYPE trueVal = max(a, b);
    std::setprecision(10);
    std::cout << "LSE         ... " << timeIt([=](){ FLOAT_TYPE z = LSE(a, b, k); }, nbTrials) << " -- Error: " << (LSE(a, b, k) - trueVal) << " (found " << LSE(a, b, k) << ")" <<  std::endl;
    std::cout << "AbsMax      ... " << timeIt([=](){ FLOAT_TYPE z = AbsMax(a, b, k); }, nbTrials) << " -- Error: " << (AbsMax(a, b, k) - trueVal) << " (found " << AbsMax(a, b, k) << ")" << std::endl;
    std::cout << "smaxLong    ... " << timeIt([=](){ FLOAT_TYPE z = smaxLong(a, b, k); }, nbTrials) << " -- Error: " << (smaxLong(a, b, k) - trueVal) << " (found " << smaxLong(a, b, k) << ")" << std::endl;
    std::cout << "smaxMini    ... " << timeIt([=](){ FLOAT_TYPE z = smaxMini(a, b, k); }, nbTrials) << " -- Error: " << (smaxMini(a, b, k) - trueVal) << " (found " << smaxMini(a, b, k) << ")" << std::endl;

    FLOAT_TYPE minBreakLSE = -1.f;
    FLOAT_TYPE minBreakAbsMax = -1.f;
    FLOAT_TYPE minBreaksmaxLong = -1.f;
    FLOAT_TYPE minBreaksmaxMini = -1.f;
    b = 0.0;
    for (int i = 0; i < 100000; i++) {
        a = FLOAT_TYPE(i);
        if (!((abs(LSE(a, b, k) - max(a, b) > 10.f) || isnan(LSE(a, b, k)) || isinf(LSE(a, b, k))))) minBreakLSE = max(a, b);
        if (!((abs(AbsMax(a, b, k) - max(a, b) > 10.f) || isnan(AbsMax(a, b, k)) || isinf(AbsMax(a, b, k))))) minBreakAbsMax = max(a, b);
        if (!((abs(smaxLong(a, b, k) - max(a, b) > 10.f) || isnan(smaxLong(a, b, k)) || isinf(smaxLong(a, b, k))))) minBreaksmaxLong = max(a, b);
        if (!((abs(smaxMini(a, b, k) - max(a, b) > 10.f) || isnan(smaxMini(a, b, k)) || isinf(smaxMini(a, b, k))))) minBreaksmaxMini = max(a, b);
    }
    std::cout << "Min break LSE: " << minBreakLSE << std::endl;
    std::cout << "Min break AbsMax: " << minBreakAbsMax << std::endl;
    std::cout << "Min break smaxLong: " << minBreaksmaxLong << std::endl;
    std::cout << "Min break smaxMini: " << minBreaksmaxMini << std::endl;

    std::cout << sizeof(float) << " " << sizeof(double) << std::endl;
}*/
