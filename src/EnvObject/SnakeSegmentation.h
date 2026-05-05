#ifndef SNAKESEGMENTATION_H
#define SNAKESEGMENTATION_H

#include "Curves/Curves.h"
#include "DataStructure/Matrix3.h"


class SnakeSegmentationParameters;
// class SnakeSegmentationImplicitParameters;
// class SnakeSegmentationExplicitParameters;

class SnakeImageField;



std::function<Vector3 (const Vector3&)> gradientFromFieldFunction(const std::function<float (const Vector3&)>& func);

class SnakeSegmentation {

public:
    SnakeSegmentation();
    SnakeSegmentation(std::shared_ptr<Curve> curve);
    SnakeSegmentation(SnakeSegmentationParameters* params, SnakeImageField* fields, std::shared_ptr<Curve> curve = std::shared_ptr<Curve>());
    SnakeSegmentation(const Curve& curve);
    SnakeSegmentation(SnakeSegmentationParameters* params, SnakeImageField* fields, const Curve& curve);
    // virtual ~SnakeSegmentation() {}
    // SnakeSegmentation(const BSpline& initialContour, const GridF& inputImage,
                      // const GridV3& inputGradient);

    std::shared_ptr<Curve> runSegmentation(int maxIterations = 100);
    std::shared_ptr<Curve> runSegmentation(std::shared_ptr<Curve> curve, int maxIterations = 100);

    Vector3 computeEnergyGradient(std::shared_ptr<Curve> contour, int index, bool usePreviousPointForInternal = true);

    Vector3 computeInternalEnergyGradient(std::shared_ptr<Curve> contour, int index, bool usePreviousPoint = true) const;
    Vector3 computeExternalEnergyGradient(std::shared_ptr<Curve> contour, int index) const;
    Vector3 computeShapeEnergyGradient(std::shared_ptr<Curve> contour, int index, bool usePreviousPoint = true) const;
    Vector3 computeGradientEnergyGradient(std::shared_ptr<Curve> contour, int index) const;

    std::shared_ptr<Curve> updateContour(std::shared_ptr<Curve> currentContour, float stepSize = 0.1f);

    float getImageAt(const Vector3& p) const;
    Vector3 getGradientImageAt(const Vector3& p) const;



    // virtual SnakeSegmentationParameters* getParameters() = 0;

// private:
    std::shared_ptr<Curve> contour;     // BSpline representing the contour
    // GridF image;         // Grayscale image grid
    // GridV3 gradientField; // Gradient field of the image

    Vector3 position = Vector3::invalid; // Initial position, attracting the whole curve

    SnakeSegmentationParameters* params;
    SnakeImageField* field;
/*
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
*/
    float stepSize = .1f;
    float currentDomainArea = 0;
    float currentIntegralOverArea = 0;
    std::vector<std::vector<float>> randomGreenCoords;

    // float internalCoefficient = 1.f;  // Coefficient for internal energy
    // float externalCoefficient = 1.f;  // Coefficient for external energy

    // float convergenceThreshold = 1e-2; // Threshold for convergence
};




class SnakeSegmentationParameters {
public:
    SnakeSegmentationParameters() {}
    ~SnakeSegmentationParameters() = default;

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


    // virtual float getImageAt(const Vector3& p) const = 0;
    // virtual Vector3 getGradientImageAt(const Vector3& p) const = 0;
};

class SnakeImageField {
public:
    SnakeImageField() {}
    virtual ~SnakeImageField() = default;

    virtual float getImage(const Vector3& p) const = 0;
    virtual Vector3 getGradient(const Vector3& p) const = 0;
};

class SnakeImageFieldExplicit : public SnakeImageField {
public:
    SnakeImageFieldExplicit() : SnakeImageFieldExplicit(GridF(),GridV3()) {}
    SnakeImageFieldExplicit(const GridF& image) : SnakeImageFieldExplicit(image, image.gradient()) {}
    SnakeImageFieldExplicit(const GridF& image, const GridV3& gradients) : SnakeImageField(), image(image), gradientField(gradients) {}

    GridF image;         // Grayscale image grid
    GridV3 gradientField; // Gradient field of the image

    virtual float getImage(const Vector3 &p) const { return image.interpolate(p); }
    virtual Vector3 getGradient(const Vector3 &p) const { return gradientField.interpolate(p); }
};

class SnakeImageFieldImplicit : public SnakeImageField {
public:
    SnakeImageFieldImplicit() : SnakeImageFieldImplicit({}, {}) {}
    SnakeImageFieldImplicit(std::function<float(const Vector3&)> imageField) : SnakeImageFieldImplicit(imageField, gradientFromFieldFunction(imageField)) {}
    SnakeImageFieldImplicit(std::function<float(const Vector3&)> imageField, std::function<Vector3(const Vector3&)> gradientField) : SnakeImageField(), imageField(imageField), gradientField(gradientField) {}

    std::function<float(const Vector3 &)> imageField; // Grayscale image grid
    std::function<Vector3(const Vector3 &)> gradientField; // Gradient field of the image

    virtual float getImage(const Vector3 &p) const { return imageField(p); }
    virtual Vector3 getGradient(const Vector3 &p) const { return gradientField(p); }
};














/*
template <class ImageField, class GradientField>
class SnakeSegmentationT {

public:
    SnakeSegmentationT();
    SnakeSegmentationT(ImageField img, GradientField grad)
        : image_(std::move(img)), grad_(std::move(grad)) {}
    virtual ~SnakeSegmentationT() {}

    BSpline runSegmentation(int maxIterations = 100);

    Vector3 computeEnergyGradient(const BSpline& contour, int index, bool usePreviousPointForInternal = true);

    Vector3 computeInternalEnergyGradient(const BSpline& contour, int index, bool usePreviousPoint = true) const;
    Vector3 computeExternalEnergyGradient(const BSpline& contour, int index) const;
    Vector3 computeShapeEnergyGradient(const BSpline& contour, int index, bool usePreviousPoint = true) const;
    Vector3 computeGradientEnergyGradient(const BSpline& contour, int index) const;

    BSpline updateContour(const BSpline& currentContour, float stepSize = 0.1f);

    float getImageAt(const Vector3& p) const {return image_(p); }
    Vector3 getGradientImageAt(const Vector3& p) const {return grad_(p); }

    // private:
    BSpline contour;     // BSpline representing the contour

    Vector3 position = Vector3::invalid(); // Initial position, attracting the whole curve

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
    std::vector<std::vector<float>> randomGreenCoords;

    ImageField image_;
    GradientField grad_;
};

struct GridImage {
    const GridF* img;
    float operator()(const Vector3& p) const noexcept {
        return img->interpolate(p);
    }
};

struct GridGrad {
    const GridV3* grad;
    Vector3 operator()(const Vector3& p) const noexcept {
        return grad->interpolate(p);
    }
};
struct ImplicitImage {
    GridF* img;
    float operator()(const Vector3& p) const noexcept { return img->at(p); }
};

struct ImplicitGrad {
    GridV3* grad;
    Vector3 operator()(const Vector3& p) const noexcept { return grad->at(p); }
};

using SnakeExplicit = SnakeSegmentationT<GridImage, GridGrad>;
using SnakeImplicit = SnakeSegmentationT<ImplicitImage, ImplicitGrad>;







template <class ImageField, class GradientField>
SnakeSegmentationT<ImageField, GradientField>::SnakeSegmentationT()
{

}

// template <class ImageField, class GradientField>
// SnakeSegmentationT<ImageField, GradientField>::SnakeSegmentationT(const BSpline &initialContour, const GridF &inputImage, const GridV3 &inputGradient)
// : SnakeSegmentationT() //, contour(initialContour), image(inputImage), gradientField(inputGradient)
// {
// gradientField = gradientField.gaussianSmooth(10.f, true, true);
// }

template <class ImageField, class GradientField>
BSpline SnakeSegmentationT<ImageField, GradientField>::runSegmentation(int maxIterations) {
    // BSpline currentContour = contour;

    float stepSize = 1.f;

    float initialTargetLength = this->targetLength;
    float initialTargetArea = this->targetArea;

    for (int iter = 0; iter < maxIterations; ++iter) {
        contour = updateContour(contour, stepSize);
        // std::cout << "Area: " << ShapeCurve(contour).computeArea() << "/" << targetArea << std::endl;

        // if (collapseFirstAndLastPoint) {
        //     currentContour.points.pop_back();
        //     currentContour.resamplePoints(currentContour.size() + 1);
        // } else {
        contour.resamplePoints();
        // }
    }
    if (collapseFirstAndLastPoint) {
        contour[-1] = contour[0];
        contour.resamplePoints(contour.size() + 1);
        contour.points.pop_back();
    }

    this->targetArea = initialTargetArea;
    this->targetLength = initialTargetLength;

    return contour;
}

template <class ImageField, class GradientField>
Vector3 SnakeSegmentationT<ImageField, GradientField>::computeEnergyGradient(const BSpline &contour, int index, bool usePreviousPointForInternal)
{
    // Compute the gradient of the total energy with respect to the control point at 'index'
    // Compute internal energy gradient
    Vector3 internalEnergyGradient = computeInternalEnergyGradient(contour, index, usePreviousPointForInternal);
    // Compute external energy gradient
    Vector3 externalEnergyGradient = computeExternalEnergyGradient(contour, index);
    // Shape gradient
    Vector3 shapeEnergyGradient = computeShapeEnergyGradient(contour, index, usePreviousPointForInternal);
    // Slope gradient
    Vector3 slopeEnergyGradient = computeGradientEnergyGradient(contour, index);

    // Combine internal and external energy gradients to get the total energy gradient
    // std::cout << "Internal: " << internalEnergyGradient.norm() << "\nExternal: " << externalEnergyGradient.norm() << "\nShape: " << shapeEnergyGradient.norm() << "\nSlope: " << slopeEnergyGradient.norm() << std::endl;
    Vector3 gradient = (internalEnergyGradient.maxMagnitude(1.f) + externalEnergyGradient.maxMagnitude(1.f) + shapeEnergyGradient.maxMagnitude(1.f) + slopeEnergyGradient.maxMagnitude(1.f)); //.normalize();
    // Vector3 gradient = internalCoefficient * internalEnergyGradient + externalCoefficient * externalEnergyGradient;

    return gradient;
}

template <class ImageField, class GradientField>
Vector3 SnakeSegmentationT<ImageField, GradientField>::computeInternalEnergyGradient(const BSpline &contour, int index, bool usePreviousPoint) const
{
    // Compute the gradient of the internal energy with respect to the control point at 'index'
    Vector3 internalEnergyGradient;
    float targetInterval = contour.length() / float(contour.size() - 1);

    Vector3 E_connectivity = Vector3();
    Vector3 E_curvature = Vector3();
    int i = index;
    int prev = i - 1;
    int next = i + 1;

    Vector3 connectVector;
    Vector3 curveVector;
    // Vector3 areaVector;

    if (i == 0 && !collapseFirstAndLastPoint) {
        E_curvature *= 0.f;
        connectVector = -(contour[next] - contour[i]);
    } else if (i == int(contour.size()) - 1 && !collapseFirstAndLastPoint) {
        E_curvature *= 0.f;
        connectVector = (contour[i] - contour[prev]);
    } else {
        connectVector = (usePreviousPoint ? (contour[i] - contour[prev]) : (contour[i] - contour[next]));
        curveVector = (contour[prev] - 2 * contour[i] + contour[next]);
        float curveLength = std::max(curveVector.length(), .1f);
        E_curvature = -2.f * curveVector / curveLength;
    }
    float connectLength = std::max(connectVector.length(), .1f);
    E_connectivity = -sign(targetInterval - connectLength) * std::pow(targetInterval - connectLength, 2) * connectVector / connectLength;

    internalEnergyGradient = curvatureCost * E_curvature + connectivityCost * E_connectivity;
    return internalEnergyGradient;
}

template <class ImageField, class GradientField>
Vector3 SnakeSegmentationT<ImageField, GradientField>::computeExternalEnergyGradient(const BSpline &contour, int index) const
{
    if (imageCost == 0) return Vector3();

    // Compute the gradient of the external energy with respect to the control point at 'index'
    Vector3 currentPoint = contour[index];

    // Get interpolated gradient from image gradient field at the current contour point
    Vector3 imageGradient = getGradientImageAt(currentPoint);

    Vector3 gradient = -imageGradient;

    if (imageInsideCoef == 0) {
        return imageCost * gradient;
    } else {
        // More complicated: compute the energy at borders and inside.
        Vector3 pos = contour[index];
        Vector3 prevPos = contour[index - 1];
        Vector3 nextPos = contour[index + 1];
        Vector3 AB = (pos - prevPos);
        Vector3 BC = (nextPos - pos);
        Vector3 areaGradientDirection = (AB.rotated90XY() + BC.rotated90XY()).normalize();
        // Area is (new upper triangle ABB' are + new lower triangle CB'B area) with Area ABC = 1/2 * |AB x() AC|
        float addedArea = -0.5f * (AB.y() * areaGradientDirection.x() - AB.x() * areaGradientDirection.y()) + 0.5f * (BC.x() * areaGradientDirection.y() - BC.y() * areaGradientDirection.x());

        std::vector<Vector3> randomPoints = ShapeCurve({pos, nextPos, prevPos}).randomPointsInside(3);
        float addedIntegral = 0;
        for (const auto& p : randomPoints) {
            addedIntegral += this->getImageAt(p);
        }

        Vector3 insideGradient = areaGradientDirection * ((currentIntegralOverArea / currentDomainArea) - ((currentIntegralOverArea + addedIntegral) / (currentDomainArea + addedArea)));
        return imageCost * (gradient * imageBordersCoef + insideGradient * imageInsideCoef) / (imageBordersCoef + imageInsideCoef);
    }
}

template <class ImageField, class GradientField>
Vector3 SnakeSegmentationT<ImageField, GradientField>::computeShapeEnergyGradient(const BSpline &contour, int index, bool usePreviousPoint) const
{
    int i = index;
    int prev = i - 1;
    int next = i + 1;

    Vector3 shapeEnergyGradient;

    if (areaCost != 0) {
        Vector3 initial = contour[i];
        ShapeCurve shape = ShapeCurve(contour);
        float area = currentDomainArea; // shape.computeArea();
        float right = shape.setPoint(i, initial + Vector3(1, 0, 0)).computeArea();
        float up = shape.setPoint(i, initial + Vector3(0, 1, 0)).computeArea();
        Vector3 areaVector = (Vector3(right - area, up - area) * sign(area - targetArea)) / float(contour.size());
        shapeEnergyGradient += areaCost * areaVector;
    }

    if (lengthCost != 0) {
        float targetInterval = targetLength / float(contour.size() - 1);

        Vector3 lengthVector;

        if (i == 0 && !collapseFirstAndLastPoint) {
            lengthVector = -(contour[next] - contour[i]);
        } else if (i == int(contour.size()) - 1 && !collapseFirstAndLastPoint) {
            lengthVector = (contour[i] - contour[prev]);
        } else {
            lengthVector = (usePreviousPoint ? (contour[i] - contour[prev]) : (contour[i] - contour[next]));
        }
        float connectLength = std::max(lengthVector.length(), .1f);
        Vector3 lengthEnergyGradient = -sign(targetInterval - connectLength) * std::pow(targetInterval - connectLength, 2) * lengthVector / connectLength;
        shapeEnergyGradient += lengthCost * lengthEnergyGradient;
    }
    return shapeEnergyGradient;
}

template <class ImageField, class GradientField>
Vector3 SnakeSegmentationT<ImageField, GradientField>::computeGradientEnergyGradient(const BSpline &contour, int index) const
{
    if (this->slopeCost != 0) {
        Vector3 gradient = getGradientImageAt(contour[index]);

        if (index == 0) {
            return slopeCost * -gradient.normalize();
        }

        float t = float(index) / float(contour.size() - 1);
        Vector3 internalEnergyGradient = slopeCost * gradient  * t;// * sign(gradient.dot(contour.getDirection(float(i) / float(contour.size() - 1))));
        return internalEnergyGradient;
    }
    return Vector3();
}

template <class ImageField, class GradientField>
BSpline SnakeSegmentationT<ImageField, GradientField>::updateContour(const BSpline &currentContour, float stepSize) {

    if (this->imageCost != 0 && this->imageInsideCoef != 0) {
        ShapeCurve contourAsRegion = ShapeCurve(currentContour);
        currentDomainArea = contourAsRegion.computeArea();
        currentIntegralOverArea = 0;
        int numberOfSamples = 100;

        if (randomGreenCoords.empty()) {
            std::vector<Vector3> randomPointsInit = contourAsRegion.randomPointsInside(numberOfSamples);
            randomGreenCoords.resize(randomPointsInit.size());
            for (size_t i = 0; i < randomGreenCoords.size(); i++) {
                randomGreenCoords[i] = computeGreenCoordinates(randomPointsInit[i], contourAsRegion);
            }
        }
        for (const auto& v : randomGreenCoords) {
            currentIntegralOverArea += getImageAt(computePointFromGreenCoordinates(v, contourAsRegion));
        }
    }


    // Initialize a new contour to be updated
    BSpline newContour = currentContour;
    int numPoints = currentContour.size();

    std::vector<Vector3> gradients(numPoints);

    float totalGradientsNorm = 0.f;

    std::vector<Vector3> internalGradients(numPoints);
    std::vector<Vector3> externalGradients(numPoints);
    std::vector<Vector3> shapeGradients(numPoints);
    std::vector<Vector3> slopeGradients(numPoints);
    float totalInternalGradients = 0.f;
    float totalExternalGradients = 0.f;
    float totalShapeGradients = 0.f;
    float totalSlopeGradients = 0.f;

    bool usePreviousPointForInternal = (random_gen::generate() > .5f ? true : false);

    for (int index = 0; index < numPoints; index++) {
        // Compute the gradient of the total energy with respect to the control point at 'index'
        // Compute internal energy gradient
        internalGradients[index] = computeInternalEnergyGradient(newContour, index, usePreviousPointForInternal).maxMagnitude(1.f);
        // Compute external energy gradient
        externalGradients[index] = computeExternalEnergyGradient(newContour, index).maxMagnitude(1.f);
        // Shape gradient
        shapeGradients[index] = computeShapeEnergyGradient(newContour, index, usePreviousPointForInternal).maxMagnitude(1.f);
        // Slope gradient
        slopeGradients[index] = computeGradientEnergyGradient(newContour, index).maxMagnitude(1.f);

        totalInternalGradients += internalGradients[index].norm();
        totalExternalGradients += externalGradients[index].norm();
        totalShapeGradients += shapeGradients[index].norm();
        totalSlopeGradients += slopeGradients[index].norm();
    }

    totalInternalGradients = (totalInternalGradients == 0 ? 1.f : totalInternalGradients / float(numPoints));
    totalExternalGradients = (totalExternalGradients == 0 ? 1.f : totalExternalGradients / float(numPoints));
    totalShapeGradients = (totalShapeGradients == 0 ? 1.f : totalShapeGradients / float(numPoints));
    totalSlopeGradients = (totalSlopeGradients == 0 ? 1.f : totalSlopeGradients / float(numPoints));

    for (int index = 0; index < numPoints; index++) {
        gradients[index] = ((internalGradients[index] / totalInternalGradients) + (externalGradients[index] / totalExternalGradients) + (shapeGradients[index] / totalShapeGradients) + (slopeGradients[index] / totalSlopeGradients));
        totalGradientsNorm += gradients[index].norm();
    }

    // std::cout << "Total energy gradient = " << totalGradientsNorm << std::endl;

    for (int i = 0; i < numPoints; ++i) {
        float normalizedStepSize = stepSize / (1.f + gradients[i].norm());
        newContour[i] -= gradients[i] * normalizedStepSize;
    }

    if (this->positionCost > 0.f && this->position.isValid()) {
        if (this->collapseFirstAndLastPoint) {
            Vector3 newCentroid = newContour.center();
            newContour.translate(this->position - newCentroid);
        } else {
            Vector3 posOnCurve = newContour.estimateClosestPos(this->position, true);
            newContour.translate(this->position - posOnCurve);
        }
    }

    auto autointersections = newContour.checkAutointersections();
    for (auto [i0, i1] : autointersections) {
        newContour[i0] = currentContour[i0];
        newContour[i0 + 1] = currentContour[i0 + 1];
        newContour[i1] = currentContour[i1];
        newContour[i1 + 1] = currentContour[i1 + 1];
    }
    return newContour;
}
*/

#endif // SNAKESEGMENTATION_H
