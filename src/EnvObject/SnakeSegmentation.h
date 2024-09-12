#ifndef SNAKESEGMENTATION_H
#define SNAKESEGMENTATION_H

#include "Utils/BSpline.h"
#include "DataStructure/Matrix3.h"

class SnakeSegmentation {

public:
    SnakeSegmentation();
    // SnakeSegmentation(const BSpline& initialContour, const GridF& inputImage,
                      // const GridV3& inputGradient);

    BSpline runSegmentation(int maxIterations = 100);

// private:
    // float computeInternalEnergy(const BSpline& contour);
    // float computeExternalEnergy(const BSpline& contour);
    // float computeEnergy(const BSpline& contour);

    Vector3 computeEnergyGradient(const BSpline& contour, int index, bool usePreviousPointForInternal = true);
    // Vector3 computeInternalEnergyGradient(const BSpline& contour, int index, bool usePreviousPoint = true);
    // Vector3 computeExternalEnergyGradient(const BSpline& contour, int index);

    Vector3 computeInternalEnergyGradient(const BSpline& contour, int index, bool usePreviousPoint = true) const;
    Vector3 computeExternalEnergyGradient(const BSpline& contour, int index) const;
    Vector3 computeShapeEnergyGradient(const BSpline& contour, int index, bool usePreviousPoint = true) const;
    Vector3 computeGradientEnergyGradient(const BSpline& contour, int index) const;

    BSpline updateContour(const BSpline& currentContour, float stepSize = 0.1f);

    // bool hasConverged(const BSpline& contour1, const BSpline& contour2);


// private:
    BSpline contour;     // BSpline representing the contour
    GridF image;         // Grayscale image grid
    GridV3 gradientField; // Gradient field of the image

    Vector3 position; // Initial position, attracting the whole curve

    float connectivityCost = 0.1f;
    float curvatureCost = 0.5f;
    float imageCost = 0.5f;
    float areaCost = 0.f;
    float lengthCost = 0.5f;
    float slopeCost = 0.f;
    float positionCost = 0.f;

    float imageBordersCoef = 1.f;
    float imageInsideCoef = 0.f;

    float targetLength = 0.f;
    float targetArea = 0.f;

    bool collapseFirstAndLastPoint = false;

    float currentDomainArea = 0;
    float currentIntegralOverArea = 0;
    std::vector<std::vector<float>> randomGreenCoords;

    // float internalCoefficient = 1.f;  // Coefficient for internal energy
    // float externalCoefficient = 1.f;  // Coefficient for external energy

    // float convergenceThreshold = 1e-2; // Threshold for convergence
};

#endif // SNAKESEGMENTATION_H
