#include "SnakeSegmentation.h"

SnakeSegmentation::SnakeSegmentation()
{

}

// SnakeSegmentation::SnakeSegmentation(const BSpline &initialContour, const GridF &inputImage, const GridV3 &inputGradient)
    // : SnakeSegmentation() //, contour(initialContour), image(inputImage), gradientField(inputGradient)
// {
    // gradientField = gradientField.gaussianSmooth(10.f, true, true);
// }

BSpline SnakeSegmentation::runSegmentation(int maxIterations) {
    BSpline currentContour = contour;

    float stepSize = .5f;

    float initialTargetLength = this->targetLength;
    float initialTargetArea = this->targetArea;

    for (int iter = 0; iter < maxIterations; ++iter) {
        float a = 0.5f + 0.5f * std::cos(float(nbCatapillars) * float(maxIterations) * 2.f * PI * float(iter) / float(maxIterations - 1));
        this->targetLength = interpolation::inv_linear(a, initialTargetLength * .15f, initialTargetLength * 2.f);
        this->targetArea = interpolation::inv_linear(a, initialTargetArea * .15f, initialTargetArea * 2.f);

        currentContour = updateContour(currentContour, stepSize);

        currentContour.resamplePoints();
    }

    contour = currentContour;

    this->targetArea = initialTargetArea;
    this->targetLength = initialTargetLength;

    return currentContour;
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
    Vector3 gradient = (internalEnergyGradient + externalEnergyGradient + shapeEnergyGradient + slopeEnergyGradient).normalize();
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

    if (i == 0 && !collapseFirstAndLastPoint) {
        E_curvature *= 0.f;
        connectVector = -(contour[next] - contour[i]);
    } else if (i == contour.size() - 1 && !collapseFirstAndLastPoint) {
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
    /*
    if (this->areaCost > 0) {
        Vector3 initial = contour[i];
        ShapeCurve shape = ShapeCurve(contour);
        float area = shape.computeArea();
        // float left = shape.setPoint(i, initial - Vector3(1, 0, 0)).computeArea();
        float right = shape.setPoint(i, initial + Vector3(1, 0, 0)).computeArea();
        // float down = shape.setPoint(i, initial - Vector3(0, 1, 0)).computeArea();
        float up = shape.setPoint(i, initial + Vector3(0, 1, 0)).computeArea();

        // areaVector = Vector3(right - left, up - down) * sign(area - targetArea);
        areaVector = Vector3(right - area, up - area).normalize() * sign(area - targetArea);
        internalEnergyGradient += areaCost * areaVector;
    }

    if (this->slopeCost != 0) {
        Vector3 gradient = gradientField.interpolate(contour[i]);
        // if (gradient.norm2() < 1e-4) {
        // internalEnergyGradient *= 0;
        // } else {
        // float diff = image.interpolate(contour[(i == 0 ? i + 1 : i)]) - image.interpolate(contour[(i == 0 ? i : i - 1)]);
        // std::cout << i << ": " << diff << "\n";
        // internalEnergyGradient += slopeCost * gradient.normalized(); // * (diff < 0 ? 0 : 1.f);
        internalEnergyGradient += slopeCost * gradient.normalize() * ((i - int(contour.size())/2) < 0 ? 1.f : -1.f);// * sign(gradient.dot(contour.getDirection(float(i) / float(contour.size() - 1))));
        // }
    }
    */
    return internalEnergyGradient;
}

Vector3 SnakeSegmentation::computeExternalEnergyGradient(const BSpline &contour, int index) const
{
    // Compute the gradient of the external energy with respect to the control point at 'index'
    Vector3 currentPoint = contour[index];

    // Get interpolated gradient from image gradient field at the current contour point
    Vector3 imageGradient = gradientField.interpolate(currentPoint);

    Vector3 gradient = -imageGradient;

    if (imageInsideCoef == 0) {
        return imageCost * gradient/*.normalize()*/;
    } else {
        // More complicated: compute the energy at borders and inside.
        Vector3 pos = contour[index];
        Vector3 prevPos = contour[index - 1];
        Vector3 nextPos = contour[index + 1];
        Vector3 AB = (pos - prevPos);
        Vector3 BC = (nextPos - pos);
        Vector3 areaGradientDirection = (AB.rotated90XY() + BC.rotated90XY()).normalize();
        // Area is (new upper triangle ABB' are + new lower triangle CB'B area) with Area ABC = 1/2 * |AB x AC|
        float addedArea = -0.5f * (AB.y * areaGradientDirection.x - AB.x * areaGradientDirection.y) + 0.5f * (BC.x * areaGradientDirection.y - BC.y * areaGradientDirection.x);

        std::vector<Vector3> randomPoints = ShapeCurve({pos, nextPos, prevPos}).randomPointsInside(3);
        float addedIntegral = 0;
        for (const auto& p : randomPoints) {
            addedIntegral += this->image.interpolate(p);
        }

        Vector3 insideGradient = areaGradientDirection * ((currentIntegralOverArea / currentDomainArea) - ((currentIntegralOverArea + addedIntegral) / (currentDomainArea + addedArea)));
        return imageCost * (gradient/*.normalize()*/ * imageBordersCoef + insideGradient/*.normalize()*/ * imageInsideCoef) / (imageBordersCoef + imageInsideCoef);
    }
}

Vector3 SnakeSegmentation::computeShapeEnergyGradient(const BSpline &contour, int index, bool usePreviousPoint) const
{
    int i = index;
    int prev = i - 1;
    int next = i + 1;

    Vector3 shapeEnergyGradient;

    if (areaCost != 0) {
        Vector3 initial = contour[i];
        ShapeCurve shape = ShapeCurve(contour);
        float area = shape.computeArea();
        float right = shape.setPoint(i, initial + Vector3(1, 0, 0)).computeArea();
        float up = shape.setPoint(i, initial + Vector3(0, 1, 0)).computeArea();
        Vector3 areaVector = Vector3(right - area, up - area).normalize() * sign(area - targetArea);
        shapeEnergyGradient += areaCost * areaVector;
    }

    if (lengthCost != 0) {
        float targetInterval = targetLength / float(contour.size() - 1);

        Vector3 lengthVector;

        if (i == 0 && !collapseFirstAndLastPoint) {
            lengthVector = -(contour[next] - contour[i]);
        } else if (i == contour.size() - 1 && !collapseFirstAndLastPoint) {
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

Vector3 SnakeSegmentation::computeGradientEnergyGradient(const BSpline &contour, int index) const
{
    if (this->slopeCost != 0) {
        Vector3 gradient = gradientField.interpolate(contour[index]);

        if (index == 0) {
            return slopeCost * -gradient.normalize();
        }

        float t = float(index) / float(contour.size() - 1);
        Vector3 currentDir = contour.getDirection(t); //(index > 0 ? (contour[index] - contour[index - 1]) : (contour[index + 1] - contour[index]));
        bool shouldGoDownward = gradient.dot(currentDir) > 0; //(index - int(contour.size())/2) < 0;
        // if (gradient.norm2() < 1e-4) {
        // internalEnergyGradient *= 0;
        // } else {
        // float diff = image.interpolate(contour[(i == 0 ? i + 1 : i)]) - image.interpolate(contour[(i == 0 ? i : i - 1)]);
        // std::cout << i << ": " << diff << "\n";
        // internalEnergyGradient += slopeCost * gradient.normalized(); // * (diff < 0 ? 0 : 1.f);
        Vector3 internalEnergyGradient = slopeCost * gradient/*.normalize()*/ /* * (shouldGoDownward ? 1.f : -1.f)*/ * t;// * sign(gradient.dot(contour.getDirection(float(i) / float(contour.size() - 1))));
        return internalEnergyGradient;
    }
    return Vector3();
}

/*
float SnakeSegmentation::computeInternalEnergy(const BSpline &contour) {
    float internalEnergy = 0.0;
    int numPoints = 100;

    float targetInterval = targetLength / float(contour.size() - 1);

    float E_connectivity = 0.f;
    float E_curvature = 0.f;
    for (int i = 0; i < contour.size() - 1; i++) {
        E_connectivity += std::pow(targetInterval- (contour[i] - contour[i + 1]).norm(), 2);
        if (i > 0) {
            E_curvature += (contour[i - 1] - 2 * contour[i] + contour[i + 1]).norm2();
        }
    }

    // float lengthError = contourLength - targetLength;
    // internalEnergy = internalEnergy / float(numPoints) + gamma * std::pow(lengthError, 2) + alpha * E_connectivity;
    internalEnergy = curvatureCost * E_curvature + lengthCost * E_connectivity;

    return internalEnergy;
}

float SnakeSegmentation::computeExternalEnergy(const BSpline &contour) {
    // Calculate external energy based on image gradients
    float externalEnergy = 0.0;
    int numPoints = contour.numPoints();

    for (int i = 0; i < numPoints; ++i) {
        Vector3 point = contour[i];
        externalEnergy += -image.interpolate(point);
    }
    // externalEnergy /= float(numPoints);
    externalEnergy = externalEnergy * imageCost;

    return externalEnergy;
}

float SnakeSegmentation::computeEnergy(const BSpline &contour)
{
    return computeInternalEnergy(contour) + computeExternalEnergy(contour);
    // return internalCoefficient * computeInternalEnergy(contour) + externalCoefficient * computeExternalEnergy(contour);
}

Vector3 SnakeSegmentation::computeEnergyGradient(const BSpline &contour, int index, bool usePreviousPointForInternal) {
    // Compute the gradient of the total energy with respect to the control point at 'index'
    // Compute internal energy gradient
    Vector3 internalEnergyGradient = computeInternalEnergyGradient(contour, index, usePreviousPointForInternal).normalize();
    // Compute external energy gradient
    Vector3 externalEnergyGradient = computeExternalEnergyGradient(contour, index).normalize();

    // Combine internal and external energy gradients to get the total energy gradient
    Vector3 gradient = (internalEnergyGradient + externalEnergyGradient).normalize();
    // Vector3 gradient = internalCoefficient * internalEnergyGradient + externalCoefficient * externalEnergyGradient;

    return gradient;
}

Vector3 SnakeSegmentation::computeInternalEnergyGradient(const BSpline &contour, int index, bool usePreviousPoint) {
    // Compute the gradient of the internal energy with respect to the control point at 'index'

    Vector3 internalEnergyGradient;
    float targetInterval = targetLength / float(contour.size() - 1);

    Vector3 E_connectivity = Vector3();
    Vector3 E_curvature = Vector3();
    int i = index;
    int prev = i - 1;
    int next = i + 1;

    Vector3 connectVector;
    Vector3 curveVector;
    Vector3 areaVector;

    if (i == 0 && !collapseFirstAndLastPoint) {
        E_curvature *= 0.f;
        connectVector = -(contour[next] - contour[i]);
    } else if (i == contour.size() - 1 && !collapseFirstAndLastPoint) {
        E_curvature *= 0.f;
        connectVector = (contour[i] - contour[prev]);
    } else {
        connectVector = (usePreviousPoint ? (contour[i] - contour[prev]) : (contour[i] - contour[next]));
        // if (connectVector.norm2() == 0) {
            // connectVector = Vector3(1e-2, 1e-2, 0) * (use)
        // }
        curveVector = (contour[prev] - 2 * contour[i] + contour[next]);
        float curveLength = std::max(curveVector.length(), .1f);
        E_curvature = -2.f * curveVector / curveLength;
    }    
    float connectLength = std::max(connectVector.length(), .1f);
    E_connectivity = -sign(targetInterval - connectLength) * std::pow(targetInterval - connectLength, 2) * connectVector / connectLength;

    internalEnergyGradient = curvatureCost * E_curvature + lengthCost * E_connectivity;

    if (this->areaCost > 0) {
        Vector3 initial = contour[i];
        ShapeCurve shape = ShapeCurve(contour);
        float area = shape.computeArea();
        // float left = shape.setPoint(i, initial - Vector3(1, 0, 0)).computeArea();
        float right = shape.setPoint(i, initial + Vector3(1, 0, 0)).computeArea();
        // float down = shape.setPoint(i, initial - Vector3(0, 1, 0)).computeArea();
        float up = shape.setPoint(i, initial + Vector3(0, 1, 0)).computeArea();

        // areaVector = Vector3(right - left, up - down) * sign(area - targetArea);
        areaVector = Vector3(right - area, up - area).normalize() * sign(area - targetArea);
        internalEnergyGradient += areaCost * areaVector;
    }

    if (this->slopeCost != 0) {
        Vector3 gradient = gradientField.interpolate(contour[i]);
        // if (gradient.norm2() < 1e-4) {
            // internalEnergyGradient *= 0;
        // } else {
        // float diff = image.interpolate(contour[(i == 0 ? i + 1 : i)]) - image.interpolate(contour[(i == 0 ? i : i - 1)]);
        // std::cout << i << ": " << diff << "\n";
        // internalEnergyGradient += slopeCost * gradient.normalized(); // * (diff < 0 ? 0 : 1.f);
        internalEnergyGradient += slopeCost * gradient.normalize() * ((i - int(contour.size())/2) < 0 ? 1.f : -1.f);// * sign(gradient.dot(contour.getDirection(float(i) / float(contour.size() - 1))));
        // }
    }
    return internalEnergyGradient;
}

Vector3 SnakeSegmentation::computeExternalEnergyGradient(const BSpline &contour, int index) {
    // Compute the gradient of the external energy with respect to the control point at 'index'

    Vector3 currentPoint = contour[index];

    // Get interpolated gradient from image gradient field at the current contour point
    Vector3 imageGradient = gradientField.interpolate(currentPoint);

    Vector3 gradient = -imageGradient;

    return imageCost * gradient;
}*/

BSpline SnakeSegmentation::updateContour(const BSpline &currentContour, float stepSize) {

    if (this->imageCost != 0 && this->imageInsideCoef != 0) {
        ShapeCurve contourAsRegion = ShapeCurve(currentContour);
        currentDomainArea = contourAsRegion.computeArea();
        currentIntegralOverArea = 0;
        int numberOfSamples = 100;

        if (randomGreenCoords.empty()) {
            std::vector<Vector3> randomPointsInit = contourAsRegion.randomPointsInside(numberOfSamples);
            randomGreenCoords.resize(randomPointsInit.size());
            for (int i = 0; i < randomGreenCoords.size(); i++) {
                randomGreenCoords[i] = computeGreenCoordinates(randomPointsInit[i], contourAsRegion);
            }
        }
        for (const auto& v : randomGreenCoords) {
            currentIntegralOverArea += image.interpolate(computePointFromGreenCoordinates(v, contourAsRegion));
        }
    }


    // Initialize a new contour to be updated
    BSpline newContour = currentContour;
    // if (positionCost != 0) {
    //     Vector3 offset = newContour.estimateClosestPos(this->position, true, 1e-1) - this->position;
    //     for (int i = 0; i < currentContour.size(); i++)
    //         newContour[i] -= offset; //.normalized();
    // }

    // Vector3 centroid = newContour.center();

    std::vector<Vector3> gradients(currentContour.size());

    float totalGradientsNorm = 0.f;

    // Iterate over each control point of the contour
    int numPoints = currentContour.size();
    for (int i = 0; i < numPoints; ++i) {
        Vector3 dir1 = computeEnergyGradient(newContour, i, (random_gen::generate() > .5f ? true : false));
        Vector3 dir2; // = computeEnergyGradient(newContour, i, false);
        Vector3 gradient = (dir1 + dir2).normalize();
        if (gradient.x != gradient.x) {
            std::cerr << "NaN found in the Snake gradient descent" << std::endl;
            gradient = computeEnergyGradient(newContour, i, true) + computeEnergyGradient(newContour, i, false);
        } else {
            totalGradientsNorm += 1.f; //gradient.norm();
            gradients[i] = gradient;
        }
    }

    if (positionCost != 0) {
        Vector3 offset = newContour.estimateClosestPos(this->position, true, 1e-1) - this->position;
        for (int i = 0; i < currentContour.size(); i++)
            gradients[i] += offset; //.normalized();
    }

    for (int i = 0; i < numPoints; ++i) {
        // if (image(newContour[i]) < 1e-5) {
            // gradients[i] = (gradients[i] + (newContour[i] - centroid).normalize() * .1f).normalize();
        // }

        // if (i == 0 || i == currentContour.size() - 1) continue;
        // if (i == 0) continue;
        newContour[i] -= gradients[i] * stepSize;
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

/*
bool SnakeSegmentation::hasConverged(const BSpline &contour1, const BSpline &contour2) {
    // Check convergence based on the movement of control points between two contours
    int numPoints = contour1.numPoints();
    float maxMovement = 0.0;

    for (int i = 0; i < numPoints; ++i) {
        // Calculate the Euclidean distance (movement) between corresponding control points
        Vector3 point1 = contour1[i];
        Vector3 point2 = contour2[i];
        float movement = (point2 - point1).norm2(); // Euclidean distance

        // Update the maximum movement observed
        if (movement > maxMovement) {
            maxMovement = movement;
        }
    }

    // Check if the maximum movement is below the convergence threshold
    std::cout << "Moved: " << maxMovement << std::endl;
    return (maxMovement < convergenceThreshold * convergenceThreshold);
}
*/
