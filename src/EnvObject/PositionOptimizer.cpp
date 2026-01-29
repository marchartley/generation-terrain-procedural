#include "PositionOptimizer.h"

#include "Graph/Pathfinding.h"
#include "Utils/Delaunay.h"

Vector3 PositionOptimizer::getHighestPosition(const Vector3 &seedPosition, [[maybe_unused]] const GridF &score, const GridV3 &gradients)
{
    return followGradient(seedPosition, gradients, 100, true);
}

Vector3 PositionOptimizer::getLowestPosition(const Vector3 &seedPosition, [[maybe_unused]] const GridF &score, const GridV3 &gradients)
{
    return followGradient(seedPosition, gradients, 100, false);
}

BSpline PositionOptimizer::trackHighestPosition(const Vector3 &seedPosition, [[maybe_unused]] const GridF &score, const GridV3 &gradients, int maxTries, bool goUp)
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
    SnakeSegmentationExplicit s;
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
}

BSpline CurveOptimizer::getExactLengthCurveFollowingGradients(const Vector3 &seedPosition, const GridF &score, const GridV3 &gradients, float targetLength)
{
    const Vector3 p = seedPosition;
    SnakeSegmentationExplicit s;
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
}

BSpline CurveOptimizer::getSkeletonCurve(const Vector3 &seedPosition, const GridF &score, const GridV3 &gradients, float targetLength)
{
    Vector3 pos = seedPosition;

    auto gradientsSmoothed = gradients.meanSmooth(5, 5, 1); //gradients.gaussianSmooth(5.f, true);
    Vector3 dir = gradientsSmoothed(pos).normalized().rotate(PI * .5f, 0, 0, 1) * targetLength * .5f;
    BSpline initialCurve = BSpline({pos - dir, pos + dir}).getPath(3);

    SnakeSegmentationExplicit s; // = SnakeSegmentationExplicit(initialCurve, score, gradients);
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

BSpline CurveOptimizer::followGradient(const Vector3 &seedPosition, [[maybe_unused]] const GridF &score, const GridV3 &gradients, int maxTries, bool goUp)
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
    SnakeSegmentationExplicit s; // = SnakeSegmentationExplicit(curve, score, gradients);
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
}

ShapeCurve AreaOptimizer::getPerimeterOptimizedShape(const Vector3 &seedPosition, const GridF &score, const GridV3 &gradients, float optmizedPerimeter)
{
    return ShapeCurve();
}

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

std::function<Vector3 (const Vector3&)> gradientFromFieldFunction(const std::function<float (const Vector3&)>& func) {
    return [=](const Vector3& pos) { float f00 = func(pos); return Vector3(func(pos + Vector3(1.f, 0.f, 0.f)) - f00, func(pos + Vector3(0.f, 1.f, 0.f)) - f00); };
}
BSpline ContinuousCurveOptimizer::getMinLengthCurveFollowingIsolevel(const Vector3 &seedPosition, const std::function<float (const Vector3 &)>& func, float minLength)
{
    const Vector3 p = seedPosition;
    SnakeSegmentationImplicit s = getSnakeForMinLengthCurveFollowingIsolevel(func, minLength);

    s.contour = BSpline({p - s.gradientField(p).rotated90XY() * minLength * .5f, p + s.gradientField(p).rotated90XY() * minLength * .5f}).resamplePoints(20);
    return s.runSegmentation(1000);
}

BSpline ContinuousCurveOptimizer::getExactLengthCurveFollowingGradients(const Vector3 &seedPosition, const std::function<float (const Vector3 &)>& func, float targetLength)
{
    const Vector3 p = seedPosition;
    SnakeSegmentationImplicit s = getSnakeForExactLengthCurveFollowingGradients(func, targetLength);
    s.contour = BSpline({p - s.gradientField(p).rotated90XY() * targetLength * .5f, p + s.gradientField(p).rotated90XY() * targetLength * .5f}).resamplePoints(20);
    return s.runSegmentation(100);
}

BSpline ContinuousCurveOptimizer::getSkeletonCurve(const Vector3 &seedPosition, const std::function<float (const Vector3 &)>& func, float targetLength)
{
    const Vector3 p = seedPosition;
    SnakeSegmentationImplicit s = getSnakeForSkeletonCurve(func, targetLength);
    Vector3 dir = s.gradientField(p).normalized().rotate(PI * .5f, 0, 0, 1) * targetLength * .5f;
    BSpline initialCurve = BSpline({p - dir, p + dir}).getPath(3);
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

BSpline ContinuousCurveOptimizer::followIsolevel(const Vector3 &seedPosition, const std::function<float (const Vector3 &)>& func, float minLength)
{
    Vector3 pos0 = seedPosition;
    // Vector3 dir0 = gradients.interpolate(pos0).normalized().cross(Vector3(0, 0, 1));
    auto gradients = gradientFromFieldFunction(func);
    Vector3 gradient;
    std::tie(pos0, gradient) = ContinuousPathOptimizer::jitterToFindPointAndGradient(pos0, Vector3::invalid(), gradients, 100, 5.f);
    if (!gradient.isValid())
        return BSpline();

    Vector3 dir0 = gradient.normalized().cross(Vector3(0, 0, 1));
    Vector3 pos1 = pos0;
    Vector3 dir1 = - dir0;

    float initialIsovalue = func(seedPosition);

    BSpline path({seedPosition});

    float totalDistance = 0.f;

    for (int i = 0; i < 100; i++) {
        // Go in the two directions at the same time

        // Start at pos0:
        std::tie(pos0, gradient) = ContinuousPathOptimizer::jitterToFindPointAndGradient(pos0, dir0, gradients, 100, 5.f);
        if (gradient.isValid()) {
            gradient.normalize();
            Vector3 newDir = gradient.cross(Vector3(0, 0, 1));
            dir0 = newDir * (dir0.dot(newDir) < 0 ? -1.f : 1.f);
            if (newDir.isValid()) {
                float newVal = func(pos0 + dir0);
                Vector3 newPos0 = ContinuousPathOptimizer::attractToIsovalue(pos0 + dir0, func, newVal, initialIsovalue, 1.f, 10);
                if (newPos0.isValid()) {
                    path.points.push_back(newPos0); // Back
                    totalDistance += (pos0 - newPos0).norm();
                    pos0 = newPos0;
                }
            }
        }

        // Then at pos1:
        std::tie(pos1, gradient) = ContinuousPathOptimizer::jitterToFindPointAndGradient(pos1, dir1, gradients, 100, 5.f);
        if (gradient.isValid()) {
            gradient.normalize();
            Vector3 newDir = gradient.cross(Vector3(0, 0, 1));
            dir1 = newDir * (dir1.dot(newDir) < 0 ? -1.f : 1.f);
            if (newDir.isValid()) {
                float newVal = func(pos1 + dir1);
                Vector3 newPos1 = ContinuousPathOptimizer::attractToIsovalue(pos1 + dir1, func, newVal, initialIsovalue, 1.f, 10);
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

BSpline ContinuousCurveOptimizer::followGradient(const Vector3 &seedPosition, const std::function<float (const Vector3 &)>& func, int maxTries, bool goUp)
{
    Vector3 pos = seedPosition;
    BSpline track;

    auto gradients = gradientFromFieldFunction(func);

    for (int i = 0; i < maxTries; i++) {
        track.points.push_back(pos);

        auto newPos = ContinuousPositionOptimizer::followGradient(pos, gradients, 1, goUp);
        if ((newPos - pos).norm2() < 1e-5) break;

        pos = newPos;
    }
    track.points.push_back(pos);
    return track;
}

SnakeSegmentationImplicit ContinuousCurveOptimizer::getSnakeForMinLengthCurveFollowingIsolevel(const std::function<float (const Vector3 &)> &func, float minLength)
{
    SnakeSegmentationImplicit s;
    s.targetLength = minLength;
    s.imageField = func;
    s.gradientField = gradientFromFieldFunction(func);

    s.areaCost = 0;
    s.collapseFirstAndLastPoint = false;
    s.curvatureCost = 0.01f;
    s.imageCost = 1.f;
    s.lengthCost = 1.f;
    s.slopeCost = 0.f; // Don't follow the slope
    return s;
}

SnakeSegmentationImplicit ContinuousCurveOptimizer::getSnakeForExactLengthCurveFollowingGradients(const std::function<float (const Vector3 &)> &func, float targetLength)
{
    SnakeSegmentationImplicit s;
    s.targetLength = targetLength;
    s.imageField = func;
    s.gradientField = gradientFromFieldFunction(func);

    s.areaCost = 0;
    s.collapseFirstAndLastPoint = false;
    s.curvatureCost = 0.1f;
    s.imageCost = 1.f;
    s.lengthCost = 10.f;
    s.slopeCost = 10.f; // Follow the slope
    return s;
}

SnakeSegmentationImplicit ContinuousCurveOptimizer::getSnakeForSkeletonCurve(const std::function<float (const Vector3 &)> &func, float targetLength)
{
    SnakeSegmentationImplicit s;
    s.imageField = func;
    s.gradientField = gradientFromFieldFunction(func);
    // s.convergenceThreshold = 1e-3;
    s.curvatureCost = 0.0f;
    s.lengthCost = 1.0f;
    s.imageCost = 1.f;
    s.targetLength = targetLength;
    return s;
}












Vector3 ContinuousPositionOptimizer::getHighestPosition(const Vector3 &seedPosition, const std::function<float (const Vector3 &)>& func)
{
    return followGradient(seedPosition, gradientFromFieldFunction(func), 100, true);
}

Vector3 ContinuousPositionOptimizer::getLowestPosition(const Vector3 &seedPosition, const std::function<float (const Vector3 &)>& func)
{
    return followGradient(seedPosition, gradientFromFieldFunction(func), 100, false);
}

BSpline ContinuousPositionOptimizer::trackHighestPosition(const Vector3 &seedPosition, const std::function<float (const Vector3 &)>& func, int maxTries, bool goUp)
{
    Vector3 pos = seedPosition;
    BSpline track;

    auto gradients = gradientFromFieldFunction(func);

    for (int i = 0; i < maxTries; i++) {
        track.points.push_back(pos);

        auto newPos = ContinuousPositionOptimizer::followGradient(pos, gradients, 1, goUp);
        if ((newPos - pos).norm2() < 1e-5) break;

        pos = newPos;
    }
    track.points.push_back(pos);
    return track;
}

Vector3 ContinuousPositionOptimizer::followGradient(const Vector3 &seedPosition, const std::function<Vector3 (const Vector3 &)>& gradients, int maxTries, bool goUp)
{
    Vector3 pos = seedPosition;
    float epsilon = 1e-8;
    float displaceFactor = 1.f * (goUp ? 1.f : -1.f);
    Vector3 previousDir;

    for (int i = 0; i < maxTries; i++) {
        if (!pos.isValid()) break;
        Vector3 grad = gradients(pos);
        if (grad.norm2() < epsilon) break;

        pos += grad.normalized() * displaceFactor;
        if (grad.dot(previousDir) < 0) {
            displaceFactor *= .5f;
        }
        previousDir = grad;
    }
    return pos;
}











std::pair<Vector3, Vector3> ContinuousPathOptimizer::jitterToFindPointAndGradient(const Vector3 &pos, const Vector3 &previousDir, const std::function<Vector3 (const Vector3 &)>& gradients, int maxTries, float jitterMaxRadius)
{
    Vector3 gradient(false);
    for (int iTry = 0; iTry < maxTries; iTry++) {
        Vector3 jitter = Vector3::random() * (jitterMaxRadius * float(iTry) / (float(maxTries)));
        if (previousDir.isValid() && jitter.dot(previousDir) < 0) continue;
        auto testPos = pos + jitter;
        gradient = gradients(testPos);
        if (gradient.isValid() && gradient.norm2() > 1e-8) {
            return {testPos, gradient};
        }
    }
    return {Vector3::invalid(), Vector3::invalid()};
}

Vector3 ContinuousPathOptimizer::attractToIsovalue(const Vector3 &pos, const std::function<float (const Vector3 &)>& func, float currentIsovalue, float targetIsovalue, float maxRectificationDistance, int nbEvaluations)
{
    float epsilon = 1e-5;
    auto gradients = gradientFromFieldFunction(func);

    if (std::abs(currentIsovalue - targetIsovalue) < epsilon) {
        Vector3 newGrad = gradients(pos).normalized();
        float bestRectificationScale = 0.f;
        float closestIso = std::numeric_limits<float>::max();
        for (int i = 0; i < nbEvaluations; i++) {
            float scale = maxRectificationDistance * float(i) / float(nbEvaluations - 1) * (currentIsovalue < targetIsovalue ? 1.f : -1.f);
            float newDiff = func(pos + newGrad * scale) - targetIsovalue;
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














ShapeCurve ContinuousAreaOptimizer::getInitialShape(const Vector3 &seedPosition, const std::function<float (const Vector3 &)>& func)
{
    ShapeCurve finalIsoline;
    Vector3 pos = seedPosition;
    ShapeCurve bestCurve;
    Vector3 jitterPos = pos;
    // Create a "curve" with maximal length as possible
    finalIsoline = ContinuousCurveOptimizer::followIsolevel(jitterPos, func, std::numeric_limits<float>::max()); //this->computeNewObjectsShapeAtPosition(jitterPos, gradients, score, directionLength).close();
    if (finalIsoline.size() > 5 && (finalIsoline.points.front() - finalIsoline.points.back()).norm2() < 3*3) {
        bestCurve = finalIsoline;
        bestCurve.closed = true;
    }
    return bestCurve;
}

ShapeCurve ContinuousAreaOptimizer::getAreaOptimizedShape(const Vector3 &seedPosition, const std::function<float (const Vector3 &)>& func, float targetArea)
{
    float fakeRadius = std::sqrt(targetArea) * .5f;

    ShapeCurve curve = ShapeCurve::circle(fakeRadius * .5f, seedPosition, 20);
    SnakeSegmentationImplicit s = getSnakeForAreaOptimizedShape(func, targetArea);
    s.contour = curve;
    BSpline result = s.runSegmentation(20);
    std::cout << result.length() << " " << ShapeCurve(result).computeArea() << " / " << s.targetArea << std::endl;
    return result;
}

SnakeSegmentationImplicit ContinuousAreaOptimizer::getSnakeForAreaOptimizedShape(const std::function<float (const Vector3 &)> &func, float targetArea)
{
    SnakeSegmentationImplicit s;
    s.imageField = func;
    s.gradientField = gradientFromFieldFunction(func);
    // s.convergenceThreshold = 1e-3;

    s.connectivityCost = 0.01f;
    s.curvatureCost = 0.0f;
    s.lengthCost = 0.0f;
    s.areaCost = 1.f;
    s.imageCost = 10.0f;
    s.targetLength = 0;
    s.targetArea = targetArea;
    s.collapseFirstAndLastPoint = true;
    s.imageInsideCoef = 1.f;
    s.imageBordersCoef = 0.f;
    return s;
}

ShapeCurve ContinuousAreaOptimizer::getPerimeterOptimizedShape(const Vector3 &seedPosition, const std::function<float (const Vector3 &)>& func, float optmizedPerimeter)
{
    return ShapeCurve();
}

SnakeSegmentationImplicit ContinuousAreaOptimizer::getSnakeForPerimeterOptimizedShape(const std::function<float (const Vector3 &)> &func, float optmizedPerimeter)
{
    return SnakeSegmentationImplicit();
}
