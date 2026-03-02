#include "SnakeSegmentation.h"

SnakeSegmentation::SnakeSegmentation()
    : SnakeSegmentation(nullptr, nullptr)
{

}

SnakeSegmentation::SnakeSegmentation(const BSpline &curve)
    : SnakeSegmentation(nullptr, nullptr, curve)
{

}

SnakeSegmentation::SnakeSegmentation(SnakeSegmentationParameters* params, SnakeImageField* fields, const BSpline& curve)
    : params(params), field(fields), contour(curve)
{

}

// SnakeSegmentation::SnakeSegmentation(const BSpline &initialContour, const GridF &inputImage, const GridV3 &inputGradient)
    // : SnakeSegmentation() //, contour(initialContour), image(inputImage), gradientField(inputGradient)
// {
    // gradientField = gradientField.gaussianSmooth(10.f, true, true);
// }

BSpline SnakeSegmentation::runSegmentation(const BSpline &curve, int maxIterations)
{
    BSpline currentContour = curve;

    float initialTargetLength = this->params->targetLength;
    float initialTargetArea = this->params->targetArea;

    for (int iter = 0; iter < maxIterations; ++iter) {

        currentContour = updateContour(currentContour, this->stepSize);
        currentContour.resamplePoints();
    }
    if (params->collapseFirstAndLastPoint && currentContour.size() > 0) {
        currentContour[-1] = currentContour[0];
        currentContour.resamplePoints(currentContour.size() + 1);
        currentContour.points.pop_back();
    }

    this->params->targetArea = initialTargetArea;
    this->params->targetLength = initialTargetLength;

    return currentContour;
}

BSpline SnakeSegmentation::runSegmentation(int maxIterations) {
    this->contour = this->runSegmentation(this->contour, maxIterations);
    return this->contour;
}


Vector3 SnakeSegmentation::computeEnergyGradient(const BSpline &contour, int index, bool usePreviousPointForInternal)
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

Vector3 SnakeSegmentation::computeInternalEnergyGradient(const BSpline &contour, int index, bool usePreviousPoint) const
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

    if (i == 0 && !params->collapseFirstAndLastPoint) {
        E_curvature *= 0.f;
        connectVector = -(contour[next] - contour[i]);
    } else if (i == int(contour.size()) - 1 && !params->collapseFirstAndLastPoint) {
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

    internalEnergyGradient = params->curvatureCost * E_curvature + params->connectivityCost * E_connectivity;
    return internalEnergyGradient;
}

Vector3 SnakeSegmentation::computeExternalEnergyGradient(const BSpline &contour, int index) const
{
    if (params->imageCost == 0) return Vector3();

    // Compute the gradient of the external energy with respect to the control point at 'index'
    Vector3 currentPoint = contour[index];

    // Get interpolated gradient from image gradient field at the current contour point
    Vector3 imageGradient = getGradientImageAt(currentPoint);

    Vector3 gradient = -imageGradient;

    if (params->imageInsideCoef == 0) {
        return params->imageCost * gradient/*.normalize()*/;
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
        if (ShapeCurve(contour).containsXY(pos + nextPos + prevPos) / 3.f) {
            addedIntegral *= -1.f;
        }

        Vector3 insideGradient = areaGradientDirection * ((currentIntegralOverArea / currentDomainArea) - ((currentIntegralOverArea + addedIntegral) / (currentDomainArea + addedArea)));
        return params->imageCost * (gradient/*.normalize()*/ * params->imageBordersCoef + insideGradient/*.normalize()*/ * params->imageInsideCoef) / (params->imageBordersCoef + params->imageInsideCoef);
    }
}

Vector3 SnakeSegmentation::computeShapeEnergyGradient(const BSpline &contour, int index, bool usePreviousPoint) const
{
    int i = index;
    int prev = i - 1;
    int next = i + 1;

    Vector3 shapeEnergyGradient;

    if (params->areaCost != 0) {
        Vector3 initial = contour[i];
        ShapeCurve shape = ShapeCurve(contour);
        float area = currentDomainArea; // shape.computeArea();
        float right = shape.setPoint(i, initial + Vector3(1, 0, 0)).computeArea();
        float up = shape.setPoint(i, initial + Vector3(0, 1, 0)).computeArea();
        Vector3 areaVector = (Vector3(right - area, up - area) * sign(area - params->targetArea)) / float(contour.size());
        shapeEnergyGradient += params->areaCost * areaVector;
    }

    if (params->lengthCost != 0) {
        float targetInterval = params->targetLength / float(contour.size() - 1);

        Vector3 lengthVector;

        if (i == 0 && !params->collapseFirstAndLastPoint) {
            lengthVector = -(contour[next] - contour[i]);
        } else if (i == int(contour.size()) - 1 && !params->collapseFirstAndLastPoint) {
            lengthVector = (contour[i] - contour[prev]);
        } else {
            lengthVector = (usePreviousPoint ? (contour[i] - contour[prev]) : (contour[i] - contour[next]));
        }
        float connectLength = std::max(lengthVector.length(), .1f);
        Vector3 lengthEnergyGradient = -sign(targetInterval - connectLength) * std::pow(targetInterval - connectLength, 2) * lengthVector / connectLength;
        shapeEnergyGradient += params->lengthCost * lengthEnergyGradient;
    }
    return shapeEnergyGradient;
}

Vector3 SnakeSegmentation::computeGradientEnergyGradient(const BSpline &contour, int index) const
{
    if (this->params->slopeCost != 0) {
        Vector3 gradient = getGradientImageAt(contour[index]);

        if (index == 0) {
            return params->slopeCost * -gradient.normalize();
        }

        float t = float(index) / float(contour.size() - 1);
        // Vector3 currentDir = contour.getDirection(t); //(index > 0 ? (contour[index] - contour[index - 1]) : (contour[index + 1] - contour[index]));
        // bool shouldGoDownward = gradient.dot(currentDir) > 0; //(index - int(contour.size())/2) < 0;
        // if (gradient.norm2() < 1e-4) {
        // internalEnergyGradient *= 0;
        // } else {
        // float diff = getImageAt(contour[(i == 0 ? i + 1 : i)]) - getImageAt(contour[(i == 0 ? i : i - 1)]);
        // std::cout << i << ": " << diff << "\n";
        // internalEnergyGradient += slopeCost * gradient.normalized(); // * (diff < 0 ? 0 : 1.f);
        Vector3 internalEnergyGradient = params->slopeCost * gradient/*.normalize()*/ /* * (shouldGoDownward ? 1.f : -1.f)*/ * t;// * sign(gradient.dot(contour.getDirection(float(i) / float(contour.size() - 1))));
        return internalEnergyGradient;
    }
    return Vector3();
}

BSpline SnakeSegmentation::updateContour(const BSpline &currentContour, float stepSize) {
    if (currentContour.empty()) return currentContour;

    ShapeCurve contourAsRegion = ShapeCurve(currentContour);
    currentDomainArea = contourAsRegion.computeArea();
    if (this->params->imageCost != 0 && this->params->imageInsideCoef != 0 && this->params->collapseFirstAndLastPoint) {
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

    #pragma omp parallel for
    for (int index = 0; index < numPoints; index++) {
        // Compute the gradient of the total energy with respect to the control point at 'index'
        // Compute internal energy gradient
        internalGradients[index] = computeInternalEnergyGradient(newContour, index, usePreviousPointForInternal).normalize(); // .maxMagnitude(1.f);
        // Compute external energy gradient
        externalGradients[index] = computeExternalEnergyGradient(newContour, index).normalize(); // .maxMagnitude(1.f);
        // Shape gradient
        shapeGradients[index] = computeShapeEnergyGradient(newContour, index, usePreviousPointForInternal).normalize(); // .maxMagnitude(1.f);
        // Slope gradient
        slopeGradients[index] = computeGradientEnergyGradient(newContour, index).normalize(); // .maxMagnitude(1.f);
    }
    for (int index = 0; index < numPoints; index++) {
        totalInternalGradients += internalGradients[index].norm();
        totalExternalGradients += externalGradients[index].norm();
        totalShapeGradients += shapeGradients[index].norm();
        totalSlopeGradients += slopeGradients[index].norm();
    }

    totalInternalGradients = (totalInternalGradients == 0 ? 1.f : totalInternalGradients / float(numPoints));
    totalExternalGradients = (totalExternalGradients == 0 ? 1.f : totalExternalGradients / float(numPoints));
    totalShapeGradients = (totalShapeGradients == 0 ? 1.f : totalShapeGradients / float(numPoints));
    totalSlopeGradients = (totalSlopeGradients == 0 ? 1.f : totalSlopeGradients / float(numPoints));

    float totalEnergyCosts = (params->curvatureCost + params->connectivityCost + params->imageCost + params->slopeCost + params->lengthCost + params->areaCost);
    if (totalEnergyCosts == 0) return currentContour;

    for (int index = 0; index < numPoints; index++) {
        gradients[index] = (
            (internalGradients[index] / totalInternalGradients) * ((params->curvatureCost + params->connectivityCost) / totalEnergyCosts) +
                            (externalGradients[index] / totalExternalGradients) * (params->imageCost / totalEnergyCosts) +
                            (shapeGradients[index] / totalShapeGradients) * ((params->lengthCost + params->areaCost) / totalEnergyCosts) +
                            (slopeGradients[index] / totalSlopeGradients) * (params->slopeCost / totalEnergyCosts)
                            );
        totalGradientsNorm += gradients[index].norm();
    }

    // std::cout << "Total energy gradient = " << totalGradientsNorm << std::endl;

    for (int i = 0; i < numPoints; ++i) {
        float normalizedStepSize = stepSize / (1.f + gradients[i].norm());
        newContour[i] -= gradients[i] * normalizedStepSize;
    }

    if (this->params->positionCost > 0.f && this->position.isValid()) {
        if (this->params->collapseFirstAndLastPoint) {
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

float SnakeSegmentation::getImageAt(const Vector3 &p) const
{
    return field->getImage(p);
}

Vector3 SnakeSegmentation::getGradientImageAt(const Vector3 &p) const
{
    return field->getGradient(p);
}





std::function<Vector3 (const Vector3&)> gradientFromFieldFunction(const std::function<float (const Vector3&)>& func) {
    return [=](const Vector3& pos) { float f00 = func(pos); return Vector3(func(pos + Vector3(1.f, 0.f, 0.f)) - f00, func(pos + Vector3(0.f, 1.f, 0.f)) - f00); };
}
