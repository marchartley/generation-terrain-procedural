#ifndef CATMULLROMSPLINE_H
#define CATMULLROMSPLINE_H

#include "Curves/Curve.h"

class CatmullRomSpline : public Curve
{
public:
    CatmullRomSpline();
    CatmullRomSpline(const std::vector<Vector3>& points);
    CatmullRomSpline(std::vector<CatmullRomSpline> subsplines);
    CatmullRomSpline(const Curve& curve);

    CLONE_FUNCTION(CatmullRomSpline)
    std::vector<Vector3> getPath(int numberOfPoints = -1) const override;
    Vector3 getPoint(float x) const override;
    Vector3 getDerivative(float x, bool normalize = false) const override;
    Vector3 getSecondDerivative(float x, bool normalize = false) const override;
    float estimateClosestTime(const Vector3& pos) const override;
    Vector3 estimateClosestPos(const Vector3& pos) const override;
    float estimateSqrDistanceFrom(const Vector3& pos) const override;
    float length() const override;

    std::vector<Vector3> getPoints() const;
    Vector3& get(int i) override;
    Vector3 get(int i) const override;

    CatmullRomSpline smooth(float factor = 1.f) const;
    CatmullRomSpline taubinSmooth(float factor = 1.f) const;

    CatmullRomSpline& setPoint(int i, const Vector3& newPos) override;

    CatmullRomSpline& resamplePoints(int newNbPoints = -1) override;

    CatmullRomSpline& reverseVertices() override;

    size_t nextID(int i);
    size_t prevID(int i);

    operator bool() const;

    std::tuple<Vector3, Vector3, Vector3> getFrenetFrame(float x) const;
    Vector3 getFrenetDirection(float x) const;
    Vector3 getFrenetNormal(float x) const;
    Vector3 getFrenetBinormal(float x) const;

    Vector3 getCenterCircle(float x) const;

    CatmullRomSpline& close() override;

    CatmullRomSpline& cleanPoints();

    Vector3 getCatmullPoint(float x) const;

    CatmullRomSpline simplifyByRamerDouglasPeucker(float epsilon, CatmullRomSpline subspline = CatmullRomSpline());

    std::pair<Vector3, Vector3> AABBox() const;

    using Curve::scale;
    CatmullRomSpline& scale(const Vector3& factor) override;
    CatmullRomSpline scaled(float factor);
    CatmullRomSpline scaled(const Vector3& factor);

    //    BSpline& grow(float increase);
    //    BSpline& shrink(float decrease);

    // CatmullRomSpline computeConvexHull() const;

    CatmullRomSpline& translate(const Vector3& translation) override;

    // std::vector<std::pair<size_t, size_t>> checkAutointersections() const;

    CatmullRomSpline& displacePointsRandomly(float maxDistance);
    CatmullRomSpline& displacePointsRandomly(const Vector3& maxDistance);
    CatmullRomSpline& displacePointsRandomlyPerlin(float maxDistance, float scale = 1.f, bool loop = false);
    CatmullRomSpline& displacePointsRandomlyPerlin(const Vector3 &maxDistance, float scale = 1.f, bool loop = false);

    virtual CatmullRomSpline& removeDuplicates() override;

    std::string toString() const override;

    std::vector<Vector3>::const_iterator begin() const override;
    std::vector<Vector3>::const_iterator end() const override;
    std::vector<Vector3>::iterator begin() override;
    std::vector<Vector3>::iterator end() override;

    std::size_t size() const;
    std::size_t numPoints() const override;
    std::size_t numVertices() const;
    bool empty() const;

    // Vector3 front() const { return (size() > 0 ? points[0] : Vector3::invalid); }
    // Vector3 back() const { return (size() > 0 ? points[size() - 1] : Vector3::invalid); }

    Vector3& operator[](size_t i);
    const Vector3& operator[](size_t i) const;

    std::string display1DPlot(int sizeX, int sizeY) const;


    Vector3 computeDerivative(float x) const;

    std::pair<Vector3, Vector3> pointAndDerivative(float x) const;
    std::tuple<Vector3, Vector3, Vector3> pointAndDerivativeAndSecondDerivative(float x) const;


    static CatmullRomSpline random(int numberOfPoints);

    void setAlpha(float newAlpha);
    float getAlpha() const;

    void addPoint(const Vector3& newPoint) override;
    CatmullRomSpline& insertPoint(int i, const Vector3& newPos) override;
    CatmullRomSpline& removePoint(int i) override;

    CatmullRomSpline& reset();


    using Curve::slice;
    std::vector<std::shared_ptr<Curve>> slice(float t) const override;


    static float CatmullNextT(const Vector3& P0, const Vector3& P1, float t_prev, float alpha);

    float alpha = 0.5f;  // alpha : 2 = very round, 1 = quite normal, 0.5 = almost linear
protected:
    std::vector<Vector3> points;
};

#endif // CATMULLROMSPLINE_H
