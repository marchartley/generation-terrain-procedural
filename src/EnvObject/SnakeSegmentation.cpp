#include "SnakeSegmentation.h"

SnakeSegmentation::SnakeSegmentation(const BSpline &initialContour, const GridF &inputImage, const GridV3 &inputGradient, float internalCoeff, float externalCoeff, float convThreshold)
    : contour(initialContour), image(inputImage), gradientField(inputGradient),
    // internalCoefficient(internalCoeff), externalCoefficient(externalCoeff),
    convergenceThreshold(convThreshold)
{
    // gradientField = gradientField.gaussianSmooth(10.f, true, true);
}

BSpline SnakeSegmentation::runSegmentation(int maxIterations) {
    BSpline currentContour = contour;
    // float realTarget = targetLength;
    // float minTargetLength = realTarget * .5f;

    float stepSize = .1f;

    int iter = 0;
    for (iter = 0; iter < maxIterations; ++iter) {
        // stepSize = random_gen::generate(2.f);
        currentContour = updateContour(currentContour, stepSize);
        /*if (maxIterations / 5 <= iter && iter <= 4 * maxIterations / 5) {
            targetLength = interpolation::inv_linear(2.f * std::abs((float(iter) / float(maxIterations - 1)) - 0.5f), minTargetLength, realTarget);
        } else {
            targetLength = realTarget;
        }*/
        // Check convergence
        /*if (hasConverged(contour, currentContour)) {
            // std::cout << "Converged at step " << iter + 1 << std::endl;
            // break;
            stepSize *= .5f;
            if (stepSize < 1e-5) {
                std::cout << "Convergence at iteration " << iter + 1 << std::endl;
                return contour;
            }
            std::cout << stepSize << std::endl;
        }*/
        currentContour.resamplePoints();
        contour = currentContour;
    }

    contour = currentContour;

    return currentContour;
}

float SnakeSegmentation::computeInternalEnergy(const BSpline &contour) {
    float internalEnergy = 0.0;
    int numPoints = 100;

    /*for (int i = 0; i < numPoints; ++i) {
        float t = float(i) / float(numPoints - 1);
        Vector3 C;
        Vector3 C_prime;
        Vector3 C_second;
        std::tie(C, C_prime, C_second) = contour.pointAndDerivativeAndSecondDerivative(t);

        // float tangentCost = C_prime.norm2();
        float curvatureCost = C_second.norm2();
        internalEnergy += beta * curvatureCost; // + alpha * tangentCost;
    }*/

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

    /*BSpline tmpCurve = contour;
    float actual = computeInternalEnergy(tmpCurve);
    float right = computeInternalEnergy(tmpCurve.setPoint(index, contour[index] + Vector3(1, 0, 0)));
    float top = computeInternalEnergy(tmpCurve.setPoint(index, contour[index] + Vector3(0, 1, 0)));

    return Vector3(right - actual, top - actual, 0);*/

    Vector3 internalEnergyGradient;
    // int numPoints = 100;

    /*for (int i = 0; i < numPoints; ++i) {
        float t = float(i) / float(numPoints - 1);
        Vector3 C;
        Vector3 C_prime;
        Vector3 C_second;
        std::tie(C, C_prime, C_second) = contour.pointAndDerivativeAndSecondDerivative(t);

        // float tangentCost = C_prime.norm2();
        float curvatureCost = C_second.norm2();
        internalEnergy += beta * curvatureCost; // + alpha * tangentCost;
    }*/

    float targetInterval = targetLength / float(contour.size() - 1);

    Vector3 E_connectivity = Vector3();
    Vector3 E_curvature = Vector3();
    int i = index;
    int prev = i - 1;
    int next = i + 1;

    Vector3 connectVector;
    Vector3 curveVector;

    if (i == 0) {
        E_curvature *= 0.f;
        connectVector = -(contour[next] - contour[i]);
    } else if (i == contour.size() - 1) {
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

    // float lengthError = contourLength - targetLength;
    // internalEnergy = internalEnergy / float(numPoints) + gamma * std::pow(lengthError, 2) + alpha * E_connectivity;
    internalEnergyGradient = curvatureCost * E_curvature + lengthCost * E_connectivity;

    // if (i == 0) {
    //     std::cout << internalEnergyGradient << std::endl;
    // }

    return internalEnergyGradient;
}

Vector3 SnakeSegmentation::computeExternalEnergyGradient(const BSpline &contour, int index) {
    // Compute the gradient of the external energy with respect to the control point at 'index'

    Vector3 currentPoint = contour[index];

    // Get interpolated gradient from image gradient field at the current contour point
    Vector3 imageGradient = gradientField.interpolate(currentPoint);

    Vector3 gradient = -imageGradient;

    return imageCost * gradient;
}

BSpline SnakeSegmentation::updateContour(const BSpline &currentContour, float stepSize) {
    // Initialize a new contour to be updated
    BSpline newContour = currentContour;

    std::vector<Vector3> gradients(currentContour.size());

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
            gradients[i] = gradient;
        }
    }

    for (int i = 0; i < numPoints; ++i) {
        newContour[i] -= gradients[i] * stepSize;
    }
    return newContour;
}

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
