#ifndef OMP_ALGO_H
#define OMP_ALGO_H

#include "DataStructure/Matrix.h"

typedef Matrix MatrixXd;
typedef Matrix VectorXd;

MatrixXd omp(const MatrixXd& D, const MatrixXd& X, const int sparsity);
std::pair<std::vector<Matrix>, std::vector<Matrix>> createDictionary(int nbSamples, int highResSize, int lowResSize);
std::pair<std::vector<Matrix>, std::vector<Matrix>> createDictionaryFromFile(std::string filename, int highResSize, int lowResSize);
Matrix flattenDictionary(const std::vector<Matrix>& dictionaries);
Matrix reconstructImage(const Matrix& coefficients, const Matrix& D, size_t imageHeight, size_t imageWidth, size_t patchSize);


#endif // OMP_ALGO_H
