#include "KarstHoleProfile.h"

KarstHoleProfile::KarstHoleProfile()
{
    this->vertices = KarstHoleProfile::createTubeProfile();
}

KarstHoleProfile::KarstHoleProfile(KarstHolePredefinedShapes shape)
{
    switch (shape) {
    case TUBE:
        this->vertices = KarstHoleProfile::createTubeProfile();
        break;
    case SOLUBLE_BED:
        this->vertices = KarstHoleProfile::createSolubleBedProfile();
        break;
    case PASSAGE:
        this->vertices = KarstHoleProfile::createPassageProfile();
        break;
    case KEYHOLE:
        this->vertices = KarstHoleProfile::createKeyholeProfile();
        break;
    case CANYON:
        this->vertices = KarstHoleProfile::createCanyonProfile();
        break;
    }
}

KarstHoleProfile::KarstHoleProfile(BSpline shape)
{
    this->vertices = shape;
}

KarstHoleProfile::KarstHoleProfile(std::vector<Vector3> shape)
{
    this->vertices = BSpline(shape);
}

KarstHoleProfile &KarstHoleProfile::rotateTowardVector(Vector3 new_dir)
{
    Vector3 forward(0, 0, 1);
    new_dir.normalize();
    float angle = std::acos(forward.dot(new_dir));
    for (Vector3& point : this->vertices.points)
        point.rotate(angle, new_dir.cross(forward).normalized());
    return *this;
}

KarstHoleProfile &KarstHoleProfile::translate(Vector3 translation)
{
    for (Vector3& point : this->vertices.points)
        point.translate(translation);
    return *this;
}

KarstHoleProfile KarstHoleProfile::interpolate(KarstHoleProfile other, BSpline path, float t)
{
    std::vector<Vector3> startingPoints = this->vertices.getPath(.1f);
    std::vector<Vector3> endingPoints = other.vertices.getPath(.1f);
    std::vector<Vector3> interpolatedPoints(startingPoints.size());
    for (size_t i = 0; i < startingPoints.size(); i++) {
        interpolatedPoints[i] = (startingPoints[i] * (1 - t) + endingPoints[i] * t);
    }
    KarstHoleProfile interpolation(interpolatedPoints);
    return interpolation.rotateTowardVector(path.getDerivative(t)).translate(path.getPoint(t));
}

KarstHoleProfile& KarstHoleProfile::setNumberOfVertices(int vertice_count)
{
    std::vector<Vector3> newPoints = this->vertices.getPath(1/(float)(vertice_count - 1));
    this->vertices = BSpline(newPoints);
    return *this;
}

KarstHoleProfile &KarstHoleProfile::setSize(float sizeX, float sizeY)
{
    /*
    Vector3 minVec(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), 0);
    Vector3 maxVec(std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), 1);
    for (const Vector3& point : this->vertices.points) {
        minVec.x = std::min(minVec.x, point.x);
        minVec.y = std::min(minVec.y, point.y);
        maxVec.x = std::max(maxVec.x, point.x);
        maxVec.y = std::max(maxVec.y, point.y);
    }
    Vector3 scaling = Vector3(sizeX, sizeY, 1) / (maxVec - minVec);*/
    float scaling = std::max(sizeX, sizeY); // Dunno what to do...
    for (Vector3& point : this->vertices.points)
        point *= scaling;
    return *this;
}









BSpline KarstHoleProfile::createTubeProfile()
{
    return BSpline({
                       {-1.,  .0, 0},
                       {-.7,  .7, 0},
                       { .0,  1., 0},
                       { .7,  .7, 0},
                       { 1.,  .0, 0},
                       { .7, -.7, 0},
                       { .0, -1., 0},
                       {-.7, -.7, 0}
                   }); //.close();
}

BSpline KarstHoleProfile::createSolubleBedProfile()
{
    return BSpline({
                       {-1.,  0., 0},
                       {-.7,  .5, 0},
                       { .7,  .5, 0},
                       { 1.,  0., 0},
                       { .7, -.5, 0},
                       {-.7, -.5, 0}
                   }); //.close();
}

BSpline KarstHoleProfile::createPassageProfile()
{
    return BSpline({
                       {-1.,  0., 0},
                       {-.5,  .5, 0},
                       { .5,  .5, 0},
                       { 1.,  0., 0},
                       { .5, -.5, 0},
                       {-.5, -.5, 0}
                   }); //.close();
}

BSpline KarstHoleProfile::createKeyholeProfile()
{
    return BSpline({
                       {-1.,  0., 0},
                       {-.5,  .5, 0},
                       { .5,  .5, 0},
                       { 1.,  0., 0},
                       { .5, -.5, 0},
                       { .2, -1., 0},
                       {-.2, -1., 0},
                       {-.5, -.5, 0}
                   }); //.close();
}

BSpline KarstHoleProfile::createCanyonProfile()
{
    return BSpline({
                       {-.5,  .5, 0},
                       {-.2,  1., 0},
                       { .2,  1., 0},
                       { .5,  .5, 0},
                       { .5, -.5, 0},
                       { .2, -1., 0},
                       {-.2, -1., 0},
                       {-.5, -.5, 0}
                   }); //.close();
}
