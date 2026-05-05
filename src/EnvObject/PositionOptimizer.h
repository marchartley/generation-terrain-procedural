#ifndef POSITIONOPTIMIZER_H
#define POSITIONOPTIMIZER_H

#include "DataStructure/Vector3.h"
#include "DataStructure/Matrix3.h"
#include "Curves/Contour.h"
#include "EnvObject/SnakeSegmentation.h"

class PathOptimizer {
public:
    static std::pair<Vector3, Vector3> jitterToFindPointAndGradient(const Vector3 &pos, const Vector3 &previousDir, const GridV3 &gradients, int maxTries, float jitterMaxRadius);
    static Vector3 attractToIsovalue(const Vector3& pos, const GridF &score, const GridV3& gradients, float currentIsovalue, float targetIsovalue, float maxRectificationDistance, int nbEvaluations);
};

class PositionOptimizer
{
public:
    static Vector3 getHighestPosition(const Vector3& seedPosition, const GridF& score, const GridV3& gradients);
    static Vector3 getLowestPosition(const Vector3& seedPosition, const GridF& score, const GridV3& gradients);

    static Polyline trackHighestPosition(const Vector3& seedPosition, const GridF& score, const GridV3& gradients, int maxTries, bool goUp);

//protected:
    static Vector3 followGradient(const Vector3& seedPosition, const GridV3& gradients, int maxTries, bool goUp);
};

class CurveOptimizer
{
public:
    static Polyline getMinLengthCurveFollowingIsolevel(const Vector3& seedPosition, const GridF& score, const GridV3& gradients, float minLength);
    static Polyline getExactLengthCurveFollowingGradients(const Vector3& seedPosition, const GridF& score, const GridV3& gradients, float targetLength);
    static Polyline getSkeletonCurve(const Vector3& seedPosition, const GridF& score, const GridV3& gradients, float targetLength);
    static Polyline followIsolevel(const Vector3& seedPosition, const GridF& score, const GridV3& gradients, float minLength);
    static Polyline followGradient(const Vector3& seedPosition, const GridF& score, const GridV3& gradients, int maxTries, bool goUp);
};

class AreaOptimizer
{
public:
    static Contour getInitialShape(const Vector3& seedPosition, const GridF& score, const GridV3& gradients);
    static Contour getAreaOptimizedShape(const Vector3& seedPosition, const GridF& score, const GridV3& gradients, float targetArea);
    static Contour getPerimeterOptimizedShape(const Vector3& seedPosition, const GridF& score, const GridV3& gradients, float optmizedPerimeter);
};



class ContinuousPathOptimizer {
public:
    static std::pair<Vector3, Vector3> jitterToFindPointAndGradient(const Vector3 &pos, const Vector3 &previousDir, const std::function<Vector3 (const Vector3&)>& gradients, int maxTries, float jitterMaxRadius);
    static Vector3 attractToIsovalue(const Vector3& pos, const std::function<float(const Vector3&)>& func, float currentIsovalue, float targetIsovalue, float maxRectificationDistance, int nbEvaluations);
};

class ContinuousPositionOptimizer
{
public:
    static Vector3 getHighestPosition(const Vector3& seedPosition, const std::function<float(const Vector3&)>& func);
    static Vector3 getLowestPosition(const Vector3& seedPosition, const std::function<float(const Vector3&)>& func);

    static Polyline trackHighestPosition(const Vector3& seedPosition, const std::function<float(const Vector3&)>& func, int maxTries, bool goUp);

    //protected:
    static Vector3 followGradient(const Vector3& seedPosition, const std::function<Vector3(const Vector3&)>& gradients, int maxTries, bool goUp);
};

class ContinuousCurveOptimizer
{
public:
    static Polyline getMinLengthCurveFollowingIsolevel(const Vector3& seedPosition, const std::function<float(const Vector3&)>& func, float minLength);
    static Polyline getExactLengthCurveFollowingGradients(const Vector3& seedPosition, const std::function<float(const Vector3&)>& func, float targetLength);
    static Polyline getSkeletonCurve(const Vector3& seedPosition, const std::function<float(const Vector3&)>& func, float targetLength);
    static Polyline followIsolevel(const Vector3& seedPosition, const std::function<float(const Vector3&)>& func, float minLength);
    static Polyline followGradient(const Vector3& seedPosition, const std::function<float(const Vector3&)>& func, int maxTries, bool goUp);

    static SnakeSegmentationParameters getSnakeForMinLengthCurveFollowingIsolevel(const std::function<float(const Vector3&)>& func, float minLength);
    static SnakeSegmentationParameters getSnakeForExactLengthCurveFollowingGradients(const std::function<float(const Vector3&)>& func, float targetLength);
    static SnakeSegmentationParameters getSnakeForSkeletonCurve(const std::function<float(const Vector3&)>& func, float targetLength);
};



class ContinuousAreaOptimizer
{
public:
    static Contour getInitialShape(const Vector3& seedPosition, const std::function<float(const Vector3&)>& func);
    static Contour getAreaOptimizedShape(const Vector3& seedPosition, const std::function<float(const Vector3&)>& func, float targetArea);
    static Contour getPerimeterOptimizedShape(const Vector3& seedPosition, const std::function<float(const Vector3&)>& func, float optmizedPerimeter);

    static SnakeSegmentationParameters getSnakeForAreaOptimizedShape(const std::function<float(const Vector3&)>& func, float targetArea);
    static SnakeSegmentationParameters getSnakeForPerimeterOptimizedShape(const std::function<float(const Vector3&)>& func, float optmizedPerimeter);
};


#endif // POSITIONOPTIMIZER_H
