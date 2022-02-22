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

KarstHoleProfile& KarstHoleProfile::setNumberOfVertices(int vertice_count)
{
    std::vector<Vector3> newPoints = this->vertices.getPath(1/(float)(vertice_count - 1));
    this->vertices = BSpline(newPoints);
    return *this;
}

KarstHoleProfile &KarstHoleProfile::setSize(float sizeX, float sizeY)
{
    Vector3 minVec(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), 0);
    Vector3 maxVec(std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), 1);
    for (const Vector3& point : this->vertices.points) {
        minVec.x = std::min(minVec.x, point.x);
        minVec.y = std::min(minVec.y, point.y);
        maxVec.x = std::max(maxVec.x, point.x);
        maxVec.y = std::max(maxVec.y, point.y);
    }
    Vector3 scaling = Vector3(sizeX, sizeY, 1) / (maxVec - minVec);
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
