#ifndef POSITIONOPTIMIZER_H
#define POSITIONOPTIMIZER_H

#include "DataStructure/Vector3.h"
#include "DataStructure/Matrix3.h"
#include "Utils/ShapeCurve.h"
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

    static BSpline trackHighestPosition(const Vector3& seedPosition, const GridF& score, const GridV3& gradients, int maxTries, bool goUp);

//protected:
    static Vector3 followGradient(const Vector3& seedPosition, const GridV3& gradients, int maxTries, bool goUp);
};

class CurveOptimizer
{
public:
    static BSpline getMinLengthCurveFollowingIsolevel(const Vector3& seedPosition, const GridF& score, const GridV3& gradients, float minLength);
    static BSpline getExactLengthCurveFollowingGradients(const Vector3& seedPosition, const GridF& score, const GridV3& gradients, float targetLength);
    static BSpline getSkeletonCurve(const Vector3& seedPosition, const GridF& score, const GridV3& gradients, float targetLength);

//protected:
    static BSpline followIsolevel(const Vector3& seedPosition, const GridF& score, const GridV3& gradients, float minLength);
    static BSpline followGradient(const Vector3& seedPosition, const GridF& score, const GridV3& gradients, int maxTries, bool goUp);
};

class AreaOptimizer
{
public:
    static ShapeCurve getInitialShape(const Vector3& seedPosition, const GridF& score, const GridV3& gradients);
    static ShapeCurve getAreaOptimizedShape(const Vector3& seedPosition, const GridF& score, const GridV3& gradients, float targetArea);
    static ShapeCurve getPerimeterOptimizedShape(const Vector3& seedPosition, const GridF& score, const GridV3& gradients, float optmizedPerimeter);
};

#endif // POSITIONOPTIMIZER_H
