#include "OMP_algo.h"

#include "Utils/Utils.h"
#include "DataStructure/Matrix3.h"
#include "DataStructure/Image.h"

using namespace std;

// Orthogonal Matching Pursuit (OMP) function
MatrixXd omp(const MatrixXd& D, const MatrixXd& X, const int sparsity) {
    size_t numAtoms = D.cols();
    size_t numSignals = X.cols();
    Matrix coefficients(numAtoms, numSignals);

    for (size_t k = 0; k < numSignals; ++k) {
        vector<float> x(X.rows());
        for (size_t i = 0; i < X.rows(); ++i) {
            x[i] = X[i][k];
        }

        vector<float> residual = x;
        vector<int> selectedAtoms;
        Matrix A(X.rows(), sparsity);

        vector<float> subCoefficients(sparsity, 0);
        for (int j = 0; j < sparsity; ++j) {
            // Compute correlations
            vector<float> correlations(numAtoms, 0);
            #pragma omp parallel for
            for (size_t i = 0; i < numAtoms; ++i) {
                if (isIn(int(i), selectedAtoms)) continue;
                #pragma omp parallel for
                for (size_t r = 0; r < residual.size(); ++r) {
                    correlations[i] += D[r][i] * residual[r];
                }
            }

            // Find the index of the atom with the maximum correlation
            int maxIndex = distance(correlations.begin(), max_element(correlations.begin(), correlations.end()));
            if (std::abs(correlations[maxIndex]) < 1e-5) break;
            selectedAtoms.push_back(maxIndex);

            // Update the matrix A with the selected atom
            for (size_t i = 0; i < A.rows(); ++i) {
                A[i][j] = D[i][maxIndex];
            }

            // Solve for the coefficients
            subCoefficients = Matrix::solve(A.leftCols(j + 1), x);

            // Update the residual
            residual = x;
            float sum = 0;
            #pragma omp parallel for collapse(2)
            for (size_t i = 0; i < A.rows(); ++i) {
                for (int l = 0; l <= j; ++l) {
                    residual[i] -= A[i][l] * subCoefficients[l];
                    sum += residual[i];
                }
            }
            // x = residual;
            // if (std::abs(sum/float(residual.size())) < 1e-5) break;
        }

        // Fill the coefficients matrix
        for (size_t i = 0; i < selectedAtoms.size(); ++i) {
            coefficients[selectedAtoms[i]][k] = subCoefficients[i];
        }
    }

    return coefficients;
}

std::pair<std::vector<Matrix>, std::vector<Matrix>> createDictionary(int nbSamples, int highResSize, int lowResSize) {
    std::vector<GridF> images(nbSamples);

    for (int i = 0; i < nbSamples; i++) {
        float angle = 2.f * PI / float(nbSamples);
        float size = highResSize;
        GridF img(size, size);
        Vector3 start = Vector3(-cos(angle), -sin(angle)) * size * .5f + Vector3(size, size) * .5f;
        Vector3 end = Vector3(cos(angle), sin(angle)) * size * .5f + Vector3(size, size) * .5f;
        img.iterateParallel([&] (const Vector3& p) {
            Vector3 startToPoint = p - start;
            Vector3 segment = end - start;
            img(p) = 1.f - startToPoint.dot(segment) / segment.norm2();
        });
        images[i] = img;
    }
    images.push_back(GridF(1, 1, 1, 1.f));

    std::vector<Matrix> bigDico(images.size());
    std::vector<Matrix> smallDico(images.size());
    #pragma omp parallel for
    for (int i = 0; i < images.size(); i++) {
        bigDico[i] = Matrix(highResSize, highResSize);
        images[i].iterateParallel([&](int x, int y, int z) {
            bigDico[i][y][x] = images[i](x, y);
        });
        smallDico[i] = Matrix(lowResSize, lowResSize);
        images[i] = images[i].resize(lowResSize, lowResSize, 1);
        images[i].iterateParallel([&](int x, int y, int z) {
            smallDico[i][y][x] = images[i](x, y);
        });
    }

    return {smallDico, bigDico};
}

std::pair<std::vector<Matrix>, std::vector<Matrix>> createDictionaryFromFile(const std::string& filename, int highResSize, int lowResSize) {
    Image initialImage;
    displayProcessTime("Reading... ", [&]() {
        initialImage = Image::readFromFile(filename);
    });
    GridF data(initialImage.colorImage.getDimensions());
    data.iterateParallel([&](size_t i) {
        data[i] = initialImage.colorImage[i].x();
    });
    std::vector<GridF> images;
    // int nbX = 200;
    // int nbY = 100;

    // int w = data.sizeX / nbX;
    // int h = data.sizeY / nbY;

    int w = highResSize, h = highResSize;
    int nbX = data.sizeX / w, nbY = data.sizeY / h;

    displayProcessTime("Splitting... ", [&]() {
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < nbX; i++) {
            for (int j = 0; j < nbY; j++) {
                if (random_gen::generate() < .5f) continue;

                GridF img = data.subset(i * w, (i+1) * w, j * h, (j+1) * h).resize(highResSize, highResSize, 1).normalized();
                // if (img.sum() / float(highResSize * highResSize) < 0.5f) continue;
                // ImageViewer::get()->addImage(img)->exec();

                #pragma omp critical
                images.push_back(img);
            }
        }
    });

    std::vector<Matrix> bigDico(images.size());
    std::vector<Matrix> smallDico(images.size());
    displayProcessTime("To dictionary (" + std::to_string(images.size()) + "... ", [&]() {
        #pragma omp parallel for
        for (int i = 0; i < images.size(); i++) {
            bigDico[i] = Matrix(highResSize, highResSize);
            images[i].iterateParallel([&](int x, int y, int z) {
                bigDico[i][y][x] = images[i](x, y);
            });
            smallDico[i] = Matrix(lowResSize, lowResSize);
            images[i] = images[i].resize(lowResSize, lowResSize, 1);
            images[i].iterateParallel([&](int x, int y, int z) {
                smallDico[i][y][x] = images[i](x, y);
            });
        }
    });
    return {bigDico, smallDico};
}


Matrix flattenDictionary(const std::vector<Matrix>& dictionaries) {
    Matrix res = Matrix(dictionaries.size(), dictionaries[0].rows() * dictionaries[0].cols());

    for (int i = 0; i < dictionaries.size(); i++) {
        std::vector<float> flat = dictionaries[i].toStdVector();
        res[i] = flat;
    }
    return res;
}

Matrix reconstructImage(const Matrix& coefficients, const Matrix& D, size_t imageHeight, size_t imageWidth, size_t patchSize) {
    Matrix reconstructed(imageHeight, imageWidth);

    auto prod = Matrix::matprod(D, coefficients).toStdVector();
    for (int i = 0; i < reconstructed.size(); i++) {
        for (int j = 0; j < reconstructed[i].size(); j++) {
            reconstructed[j][i] = prod[j * imageHeight + i];
        }
    }
    return reconstructed;
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
