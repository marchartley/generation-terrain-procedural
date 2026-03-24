#ifndef OMP_ALGO_H
#define OMP_ALGO_H

#include "DataStructure/Matrix3.h"


struct OMPResult
{
    std::vector<int> support;     // selected atom indices
    std::vector<double> coeffs;    // coefficients for selected atoms
    GridF approximation;          // reconstructed signal
    GridF residual;               // final residual
};

class OMP {
public:
    OMP();

    enum class DictionaryType { Dirac, DCT, Gaussian };

    GridF getLargeReconstruction(const GridF& input, const Vector3i& finalDimensions = Vector3i::invalid);
    void createDictionaryFromImages(const std::vector<std::string>& paths);
    // void createDictionary(int width, int height, DictionaryType type);

    GridF reconstructImage(const OMPResult& result, const std::vector<GridF>& dictionary);

    int tileSize = 10;
    int upscaleFactor = 2;
    int nbPrimitives = 5;
    int overlap = 2;
    int maxAtoms = 2000;

protected:
    std::vector<double> solveLinearSystem(std::vector<std::vector<double>> M, std::vector<double> b);
    GridF combineAtoms(const std::vector<GridF>& dictionary, const std::vector<int>& support, const std::vector<double>& coeffs);
    std::vector<double> solveLeastSquaresActiveSet(const GridF& y, const std::vector<GridF>& dictionary, const std::vector<int>& support);
    OMPResult orthogonalMatchingPursuit(const GridF& y, const std::vector<GridF>& dictionary, int maxAtoms, double residualTolerance = 1e-6f, bool debug = false);

    std::vector<GridF> createDictionaryFromImage(int width, int height, const GridF& image);
    /*
    void createDictionary_Dirac(int width, int height);
    void createDictionary_DCT(int width, int height);
    void createDictionary_Gaussian(int width, int height, const std::vector<double>& sigmas);
    */
    void normalizeDictionary(std::vector<GridF>& dict);

    std::vector<GridF> smallDictionary;
    std::vector<GridF> bigDictionary;
};
#endif // OMP_ALGO_H
