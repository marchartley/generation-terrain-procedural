#ifndef SNAKESEGMENTATIONIMPLICIT_H
#define SNAKESEGMENTATIONIMPLICIT_H

#include "Utils/BSpline.h"
#include "Utils/ShapeCurve.h"
#include "EnvObject/SnakeSegmentation.h"
#include <functional>

class SnakeSegmentationImplicit : public SnakeSegmentation
{
public:
    SnakeSegmentationImplicit();


    float getImageAt(const Vector3& p) const;
    Vector3 getGradientImageAt(const Vector3& p) const;
    /*
    Vector3 computeEnergyGradient(const BSpline& contour, int index, bool usePreviousPointForInternal = true);
    Vector3 computeInternalEnergyGradient(const BSpline& contour, int index, bool usePreviousPoint = true) const;
    Vector3 computeExternalEnergyGradient(const BSpline& contour, int index) const;
    Vector3 computeShapeEnergyGradient(const BSpline& contour, int index, bool usePreviousPoint = true) const;
    Vector3 computeGradientEnergyGradient(const BSpline& contour, int index) const;

    BSpline updateContour(const BSpline& currentContour, float stepSize = 0.1f);

    // private:
    // BSpline contour;     // BSpline representing the contour
    std::function<float(const Vector3 &)> imageField; // Grayscale image grid
    std::function<Vector3(const Vector3 &)> gradientField; // Gradient field of the image

    // Vector3 position; // Initial position, attracting the whole curve

    float connectivityCost = 0.0f;
    float curvatureCost = 0.0f;
    float imageCost = 0.0f;
    float areaCost = 0.f;
    float lengthCost = 0.0f;
    float slopeCost = 0.f;

    float positionCost = 0.f;
    int nbCatapillars = 0;

    float imageBordersCoef = 1.f;
    float imageInsideCoef = 0.f;

    float targetLength = 0.f;
    float targetArea = 0.f;

    bool collapseFirstAndLastPoint = false;

    float currentDomainArea = 0;
    float currentIntegralOverArea = 0;
    std::vector<std::vector<float>> randomGreenCoords;*/

    // float internalCoefficient = 1.f;  // Coefficient for internal energy
    // float externalCoefficient = 1.f;  // Coefficient for external energy

    // float convergenceThreshold = 1e-2; // Threshold for convergence

    std::function<float(const Vector3 &)> imageField; // Grayscale image grid
    std::function<Vector3(const Vector3 &)> gradientField; // Gradient field of the image
};

#endif // SNAKESEGMENTATIONIMPLICIT_H
