#include "PositionOptimizer.h"

#include "Graph/Pathfinding.h"
#include "Utils/Delaunay.h"

Vector3 PositionOptimizer::getHighestPosition(const Vector3 &seedPosition, const GridF &score, const GridV3 &gradients)
{
    return followGradient(seedPosition, gradients, 100, true);
}

Vector3 PositionOptimizer::getLowestPosition(const Vector3 &seedPosition, const GridF &score, const GridV3 &gradients)
{
    return followGradient(seedPosition, gradients, 100, false);
}

BSpline PositionOptimizer::trackHighestPosition(const Vector3 &seedPosition, const GridF &score, const GridV3 &gradients, int maxTries, bool goUp)
{
    Vector3 pos = seedPosition;
    BSpline track;

    for (int i = 0; i < maxTries; i++) {
        track.points.push_back(pos);

        auto newPos = PositionOptimizer::followGradient(pos, gradients, 1, goUp);
        if ((newPos - pos).norm2() < 1e-5) break;

        pos = newPos;
    }
    track.points.push_back(pos);
    return track;
}

Vector3 PositionOptimizer::followGradient(const Vector3 &seedPosition, const GridV3 &gradients, int maxTries, bool goUp)
{
    Vector3 pos = seedPosition;
    float epsilon = 1e-8;
    float displaceFactor = 1.f * (goUp ? 1.f : -1.f);
    Vector3 previousDir;

    for (int i = 0; i < maxTries; i++) {
        if (!pos.isValid()) break;
        Vector3 grad = gradients.interpolate(pos);
        if (grad.norm2() < epsilon) break;

        pos += grad.normalized() * displaceFactor;
        if (grad.dot(previousDir) < 0) {
            displaceFactor *= .5f;
        }
        previousDir = grad;
    }
    return pos;
}

BSpline CurveOptimizer::getMinLengthCurveFollowingIsolevel(const Vector3 &seedPosition, const GridF &score, const GridV3 &gradients, float minLength)
{
    const Vector3 p = seedPosition;
    SnakeSegmentation s;
    s.contour = BSpline({p - gradients(p).rotated90XY() * minLength * .5f, p + gradients(p).rotated90XY() * minLength * .5f}).resamplePoints(20);
    s.targetLength = minLength;
    s.gradientField = gradients;
    s.image = score;

    s.areaCost = 0;
    s.collapseFirstAndLastPoint = false;
    s.curvatureCost = 0.01f;
    s.imageCost = 1.f;
    s.lengthCost = 1.f;
    s.slopeCost = 0.f; // Don't follow the slope

    return s.runSegmentation(1000);

    /*
    int maxTries = 5;

    Vector3 pos = seedPosition;

    BSpline curve;


    int bestDiffDownhill = maxTries;
    Vector3 bestPosDownhill;
    int bestDiffUphill = maxTries;
    Vector3 bestPosUphill;

    // Start by going downhill, there are much more chances that we find a best solution... (I guess...)
    for (int i = 0; i < maxTries; i++) {
        curve = CurveOptimizer::followIsolevel(pos, score, gradients, minLength);
        // Check if we fullfill minLength condition
        float length = curve.length();
        float diff = length - minLength;
        if (diff < 0 || length == 0) {
            Vector3 gradient;
            std::tie(pos, gradient) = PathOptimizer::jitterToFindPointAndGradient(pos, Vector3(false), gradients, 5, 2.f);
            if (pos.isValid())
                pos += gradient.normalized() * -1.f;
            else
                break;

        } else {
            bestDiffDownhill = i;
            bestPosDownhill = pos;
            break;
        }
    }
    // Now try uphill
    pos = seedPosition;
    for (int i = 0; i < bestDiffDownhill; i++) {
        curve = CurveOptimizer::followIsolevel(pos, score, gradients, minLength);
        // Check if we fullfill minLength condition
        float length = curve.length();
        float diff = length - minLength;
        if (diff < 0 || length == 0) {
            Vector3 gradient;
            std::tie(pos, gradient) = PathOptimizer::jitterToFindPointAndGradient(pos, Vector3(false), gradients, 5, 2.f);
            if (pos.isValid())
                pos += gradient.normalized() * 1.f;
            else
                break;
        } else {
            bestDiffUphill = i;
            bestPosUphill = pos;
            break;
        }
    }

    curve = CurveOptimizer::followIsolevel((bestDiffDownhill < bestDiffUphill ? bestPosDownhill : bestPosUphill), score, gradients, minLength);
    return curve;*/
}

BSpline CurveOptimizer::getExactLengthCurveFollowingGradients(const Vector3 &seedPosition, const GridF &score, const GridV3 &gradients, float targetLength)
{
    const Vector3 p = seedPosition;
    SnakeSegmentation s;
    s.contour = BSpline({p - gradients(p).rotated90XY() * targetLength * .5f, p + gradients(p).rotated90XY() * targetLength * .5f}).resamplePoints(20);
    s.targetLength = targetLength;
    s.gradientField = gradients;
    s.image = score;

    s.areaCost = 0;
    s.collapseFirstAndLastPoint = false;
    s.curvatureCost = 0.1f;
    s.imageCost = 1.f;
    s.lengthCost = 10.f;
    s.slopeCost = 10.f; // Follow the slope

    return s.runSegmentation(1000);
    /*
    int maxTries = 2.f * targetLength;
    int tries = 0;

    Vector3 pos = seedPosition;

    BSpline curve;

    auto curveA = CurveOptimizer::followGradient(pos, score, gradients, 5, true);
    auto curveB = CurveOptimizer::followGradient(pos, score, gradients, 5, false);
    curve = curveB.reverseVertices();
    for (const auto& p : curveA) {
        curve.points.push_back(p);
    }
    curve.removeDuplicates();
    curve.resamplePoints(std::ceil(targetLength));
    float length = curve.length();

    while (std::abs(length - targetLength) > 1.f && tries < maxTries) {
        float diff = length - targetLength;
        for (int iPath = 0; iPath < curve.size(); iPath++) {
            auto p = curve[iPath];
            float t = float(iPath) / float(curve.size() - 1);
            float scaleOnGradient = (t - .5f) * 2.f * std::abs(diff) / float(curve.size());
            Vector3 gradient;
            std::tie(p, gradient) = PathOptimizer::jitterToFindPointAndGradient(p, Vector3(false), gradients, 5, 2.f);
            if (p.isValid() && gradient.isValid())
                p += gradient.normalized() * scaleOnGradient * sign(diff);
            else
                continue;
            curve[iPath] = p;
        }
        curve.resamplePoints();
        length = curve.length();
        tries++;
    }
    return curve.resamplePoints();*/
}

BSpline CurveOptimizer::getSkeletonCurve(const Vector3 &seedPosition, const GridF &score, const GridV3 &gradients, float targetLength)
{
    Vector3 pos = seedPosition;

    auto gradientsSmoothed = gradients.gaussianSmooth(5.f, true);
    Vector3 dir = gradientsSmoothed(pos).normalized().rotate(PI * .5f, 0, 0, 1) * targetLength * .5f;
    BSpline initialCurve = BSpline({pos - dir, pos + dir}).getPath(3);

    SnakeSegmentation s; // = SnakeSegmentation(initialCurve, score, gradients);
    s.contour = initialCurve;
    s.image = score;
    s.gradientField = gradients;
    // s.convergenceThreshold = 1e-3;
    s.curvatureCost = 0.0f;
    s.lengthCost = 1.0f;
    s.imageCost = 1.f;
    s.targetLength = targetLength;
    s.contour = initialCurve;
    int nbIterations = 10;
    for (int i = 0; i < nbIterations; i++) {
        int nbCatapillars = 3;
        float a = 0.5f + 0.5f * std::cos(float(nbCatapillars) * float(nbIterations) * 2.f * PI * float(i) / float(nbIterations - 1));
        s.targetLength = interpolation::inv_linear(a, targetLength * .5f, targetLength);
        initialCurve = s.runSegmentation(40);
        s.contour.resamplePoints(s.contour.size() + 1);
    }
    initialCurve = s.runSegmentation(50);
    if (initialCurve.size() > 1 && (initialCurve[0] - initialCurve[-1]).norm2() < std::pow(targetLength * .3f, 2)) {
        initialCurve = BSpline();
    }
    return initialCurve;
}

BSpline CurveOptimizer::followIsolevel(const Vector3 &seedPosition, const GridF &score, const GridV3 &gradients, float minLength)
{
    Vector3 pos0 = seedPosition;
    // Vector3 dir0 = gradients.interpolate(pos0).normalized().cross(Vector3(0, 0, 1));
    Vector3 gradient;
    std::tie(pos0, gradient) = PathOptimizer::jitterToFindPointAndGradient(pos0, Vector3::invalid(), gradients, 100, 5.f);
    if (!gradient.isValid())
        return BSpline();

    Vector3 dir0 = gradient.normalized().cross(Vector3(0, 0, 1));
    Vector3 pos1 = pos0;
    Vector3 dir1 = - dir0;

    float initialIsovalue = score.interpolate(seedPosition);

    BSpline path({seedPosition});

    float totalDistance = 0.f;

    for (int i = 0; i < 100; i++) {
        // Go in the two directions at the same time

        // Start at pos0:
        std::tie(pos0, gradient) = PathOptimizer::jitterToFindPointAndGradient(pos0, dir0, gradients, 100, 5.f);
        if (gradient.isValid()) {
            gradient.normalize();
            Vector3 newDir = gradient.cross(Vector3(0, 0, 1));
            dir0 = newDir * (dir0.dot(newDir) < 0 ? -1.f : 1.f);
            if (newDir.isValid()) {
                float newVal = score.interpolate(pos0 + dir0);
                Vector3 newPos0 = PathOptimizer::attractToIsovalue(pos0 + dir0, score, gradients, newVal, initialIsovalue, 1.f, 10);
                if (newPos0.isValid()) {
                    path.points.push_back(newPos0); // Back
                    totalDistance += (pos0 - newPos0).norm();
                    pos0 = newPos0;
                }
            }
        }

        // Then at pos1:
        std::tie(pos1, gradient) = PathOptimizer::jitterToFindPointAndGradient(pos1, dir1, gradients, 100, 5.f);
        if (gradient.isValid()) {
            gradient.normalize();
            Vector3 newDir = gradient.cross(Vector3(0, 0, 1));
            dir1 = newDir * (dir1.dot(newDir) < 0 ? -1.f : 1.f);
            if (newDir.isValid()) {
                float newVal = score.interpolate(pos1 + dir1);
                Vector3 newPos1 = PathOptimizer::attractToIsovalue(pos1 + dir1, score, gradients, newVal, initialIsovalue, 1.f, 10);
                if (newPos1.isValid()) {
                    path.points.insert(path.points.begin(), newPos1); // Front
                    totalDistance += (pos1 - newPos1).norm();
                    pos1 = newPos1;
                }
            }
        }

        // Check for a loop: P0 == Pn || P0 == Pn-1 || P1 == Pn
        float loopEpsilon = 5.f;
        if (path.size() > 5 && (path.points.front() - path.points.back()).norm2() < loopEpsilon){
            break; // Got back close to beginning
        }
        if (minLength < totalDistance || path.size() > 5000) {
            break;
        }
    }
    // Run check for detecting spiral (This is my own algo, I don't think it's a good one...)
    AABBox boundingBox(path.points);
    float ratioAreaPerimeter = (std::pow(boundingBox.dimensions().maxComp(), 2) / totalDistance);
    float ratioLimit = .5f * totalDistance / (2.f * PI); // Approximatively the ratio for a circle...
    if (ratioAreaPerimeter < ratioLimit) return BSpline();
    return path;
}

BSpline CurveOptimizer::followGradient(const Vector3 &seedPosition, const GridF &score, const GridV3 &gradients, int maxTries, bool goUp)
{
    Vector3 pos = seedPosition;
    BSpline track;

    for (int i = 0; i < maxTries; i++) {
        track.points.push_back(pos);

        auto newPos = PositionOptimizer::followGradient(pos, gradients, 1, goUp);
        if ((newPos - pos).norm2() < 1e-5) break;

        pos = newPos;
    }
    track.points.push_back(pos);
    return track;
}

ShapeCurve AreaOptimizer::getInitialShape(const Vector3 &seedPosition, const GridF &score, const GridV3 &gradients)
{
    ShapeCurve finalIsoline;
    Vector3 pos = seedPosition;
    ShapeCurve bestCurve;
    Vector3 jitterPos = pos;
    // Create a "curve" with maximal length as possible
    finalIsoline = CurveOptimizer::followIsolevel(jitterPos, score, gradients, std::numeric_limits<float>::max()); //this->computeNewObjectsShapeAtPosition(jitterPos, gradients, score, directionLength).close();
    if (finalIsoline.size() > 5 && (finalIsoline.points.front() - finalIsoline.points.back()).norm2() < 3*3) {
        bestCurve = finalIsoline;
        bestCurve.closed = true;
    }
    return bestCurve;
}

ShapeCurve AreaOptimizer::getAreaOptimizedShape(const Vector3 &seedPosition, const GridF &score, const GridV3 &gradients, float targetArea)
{
    float fakeRadius = std::sqrt(targetArea) * .5f;
    float fakeArea = PI * fakeRadius * fakeRadius;

    ShapeCurve curve = ShapeCurve::circle(fakeRadius * .5f, seedPosition, 20);
    SnakeSegmentation s; // = SnakeSegmentation(curve, score, gradients);
    s.contour = curve;
    s.image = score;
    s.gradientField = gradients;
    // s.convergenceThreshold = 1e-3;

    s.connectivityCost = 0.01f;
    s.curvatureCost = 0.0f;
    s.lengthCost = 0.0f;
    s.areaCost = 1.f;
    s.imageCost = 10.0f;
    s.targetLength = 0;
    s.targetArea = fakeArea;
    s.collapseFirstAndLastPoint = true;
    s.imageInsideCoef = 1.f;
    s.imageBordersCoef = 0.f;

    BSpline result = s.runSegmentation(200);
    std::cout << result.length() << " " << ShapeCurve(result).computeArea() << " / " << s.targetArea << std::endl;
    return result;

    /*


    Vector3 currentSeedPos = seedPosition;
    float maxError = 5.f;
    int maxTries = 100;
    ShapeCurve finalCurve;
    float moveFactor = 1.f;
    bool currentlyAreaGettingSmaller = true;


    Vector3 pos = seedPosition;

    while (maxTries > 0) {
        ShapeCurve curve = AreaOptimizer::getInitialShape(currentSeedPos, score, gradients).resamplePoints(50);
        SnakeSegmentation s; // = SnakeSegmentation(curve, score, gradients);
        s.contour = curve;
        s.image = score;
        s.gradientField = gradients;
        // s.convergenceThreshold = 1e-3;
        s.curvatureCost = 0.01f;
        s.lengthCost = 0.0f;
        s.areaCost = 10.f;
        s.imageCost = .0f;
        s.targetLength = 0;
        float fakeRadius = std::sqrt(targetArea) * .5f;
        float fakeArea = PI * fakeRadius * fakeRadius;
        s.targetArea = fakeArea;
        s.contour = curve;
        s.collapseFirstAndLastPoint = true;

        curve = s.runSegmentation(100);

        if (curve.size() > 0) {
            // std::cout << "Shape area: " << curve.computeArea() << " - Target: " << s.targetArea << " (" << fakeRadius << "x" << fakeRadius << ")" << std::endl;
            return curve;
        }

        Vector3 gradient = gradients.interpolate(currentSeedPos).normalized();
        if (!gradient.isValid()) break;
        currentSeedPos = currentSeedPos + gradient * 2.f;
        maxTries--;

    }
    return ShapeCurve();
    */
/*
    // We will move only in the direction of the gradient, since we want to optimize the isolevel.
    // And we know that higher isolevel => lower area while lower isolevel => higher area.
    // So isolevel gradient proportional to -area gradient.
    while (maxTries > 0) {
        ShapeCurve curve = AreaOptimizer::getInitialShape(currentSeedPos, score, gradients);
        Vector3 gradient = gradients.interpolate(currentSeedPos).normalized();
        if (!gradient.isValid()) break;

        if (curve.size() == 0) {
            // The isocontour is too big, we didn't manage to do a full circle.
            currentSeedPos = currentSeedPos + gradient * 2.f;
        } else {
            float area = curve.computeArea();

            float diff = targetArea - area; // < 0 means curve too big, > 0 means curve too small
            finalCurve = curve;
            if (std::abs(diff) < maxError) break;
            currentSeedPos = currentSeedPos + gradient * (diff > 0 ? -1.f : 1.f) * moveFactor;

            if (currentlyAreaGettingSmaller != (diff > 0)) {
                currentlyAreaGettingSmaller = !currentlyAreaGettingSmaller;
                moveFactor *= .5f;
            }
        }
        maxTries--;
    }
    return finalCurve;
*/
}

ShapeCurve AreaOptimizer::getPerimeterOptimizedShape(const Vector3 &seedPosition, const GridF &score, const GridV3 &gradients, float optmizedPerimeter){}

std::pair<Vector3, Vector3> PathOptimizer::jitterToFindPointAndGradient(const Vector3& pos, const Vector3& previousDir, const GridV3& gradients, int maxTries, float jitterMaxRadius)
{
    Vector3 gradient(false);
    for (int iTry = 0; iTry < maxTries; iTry++) {
        Vector3 jitter = Vector3::random() * (jitterMaxRadius * float(iTry) / (float(maxTries)));
        if (previousDir.isValid() && jitter.dot(previousDir) < 0) continue;
        auto testPos = pos + jitter;
        gradient = gradients.interpolate(testPos);
        if (gradient.isValid() && gradient.norm2() > 1e-8) {
            return {testPos, gradient};
        }
    }
    return {Vector3::invalid(), Vector3::invalid()};
}

Vector3 PathOptimizer::attractToIsovalue(const Vector3 &pos, const GridF& score, const GridV3 &gradients, float currentIsovalue, float targetIsovalue, float maxRectificationDistance, int nbEvaluations)
{
    float epsilon = 1e-5;
    if (std::abs(currentIsovalue - targetIsovalue) < epsilon) {
        Vector3 newGrad = gradients.interpolate(pos).normalized();
        float bestRectificationScale = 0.f;
        float closestIso = std::numeric_limits<float>::max();
        for (int i = 0; i < nbEvaluations; i++) {
            float scale = maxRectificationDistance * float(i) / float(nbEvaluations - 1) * (currentIsovalue < targetIsovalue ? 1.f : -1.f);
            float newDiff = score.interpolate(pos + newGrad * scale) - targetIsovalue;
            if (std::abs(newDiff) < closestIso) {
                closestIso = std::abs(newDiff);
                bestRectificationScale = scale;
                if (closestIso < epsilon)
                    break;
            }
        }
        return pos + newGrad * bestRectificationScale;
    }
    return pos;
}
