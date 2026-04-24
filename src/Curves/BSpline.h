#ifndef BSPLINE_H
#define BSPLINE_H

#include "Curves/Curve.h"

class BSpline : public Curve
{
public:
    BSpline();
    BSpline(const std::vector<Vector3>& points);
    BSpline(std::vector<BSpline> subsplines);

    std::vector<Vector3> getPath(int numberOfPoints = -1) const;
    Vector3 getPoint(float x) const;
    Vector3 getPoint(float x, const Vector3& a, const Vector3& b) const;
    Vector3 getDerivative(float x, bool normalize = false) const;
    Vector3 getSecondDerivative(float x, bool normalize = false) const;
    float estimateClosestTime(const Vector3& pos, float epsilon = 1e-4, float nbChecksFactor = 2.f, float earlyExitThreshold = 1e-3) const;
    Vector3 estimateClosestPos(const Vector3& pos, bool useNativeShape = false, float epsilon = 1e-3) const;
    float estimateSqrDistanceFrom(const Vector3& pos, bool useNativeShape = false, float epsilon = 1e-3) const;
    float length() const;

    std::vector<Vector3> getPoints() const { return this->points; }

    BSpline smooth(float factor = 1.f) const;
    BSpline taubinSmooth(float factor = 1.f) const;

    BSpline& setPoint(int i, const Vector3& newPos);

    BSpline& resamplePoints(int newNbPoints = -1);

    BSpline& reverseVertices();

    size_t nextID(int i) { return (i + 1 + this->points.size()) % this->points.size(); }
    size_t prevID(int i) { return (i - 1 + this->points.size()) % this->points.size(); }

    operator bool() const { return (this->points.size() > 0); }

    std::tuple<Vector3, Vector3, Vector3> getFrenetFrame(float x) const;
    Vector3 getFrenetDirection(float x) const;
    Vector3 getFrenetNormal(float x) const;
    Vector3 getFrenetBinormal(float x) const;

    Vector3 getCenterCircle(float x) const;

    Vector3 center() const;

    BSpline& close();

    BSpline& cleanPoints();

    Vector3 getCatmullPoint(float x) const;

    BSpline simplifyByRamerDouglasPeucker(float epsilon, BSpline subspline = BSpline());

    std::pair<Vector3, Vector3> AABBox() const;

    BSpline& scale(float factor);
    BSpline& scale(const Vector3& factor);
    BSpline scaled(float factor);
    BSpline scaled(const Vector3& factor);

    //    BSpline& grow(float increase);
    //    BSpline& shrink(float decrease);

    BSpline computeConvexHull() const;

    BSpline& translate(const Vector3& translation);

    std::vector<std::pair<size_t, size_t>> checkAutointersections() const;

    BSpline& displacePointsRandomly(float maxDistance);
    BSpline& displacePointsRandomly(const Vector3& maxDistance);
    BSpline& displacePointsRandomlyPerlin(float maxDistance, float scale = 1.f, bool loop = false);
    BSpline& displacePointsRandomlyPerlin(const Vector3 &maxDistance, float scale = 1.f, bool loop = false);

    virtual BSpline& removeDuplicates();

    std::string toString() const;

    auto begin() const { return points.begin(); }
    auto end() const { return points.end(); }
    auto begin() { return points.begin(); }
    auto end() { return points.end(); }
    std::size_t size() const { return end() - begin(); }
    std::size_t numPoints() const { return size(); }
    std::size_t numVertices() const { return size(); }
    bool empty() const { return begin() == end(); }

    Vector3 front() const { return (size() > 0 ? points[0] : Vector3::invalid); }
    Vector3 back() const { return (size() > 0 ? points[size() - 1] : Vector3::invalid); }

    Vector3& operator[](size_t i);
    const Vector3& operator[](size_t i) const;

    std::string display1DPlot(int sizeX, int sizeY) const;


    Vector3 computeDerivative(float x) const;

    std::pair<Vector3, Vector3> pointAndDerivative(float x) const;
    std::tuple<Vector3, Vector3, Vector3> pointAndDerivativeAndSecondDerivative(float x) const;


    static BSpline random(int numberOfPoints);

    void setAlpha(float newAlpha) { this->alpha = newAlpha; }
    float getAlpha() const { return alpha; }

    void addPoint(const Vector3& newPoint);
    BSpline& insertPoint(int i, const Vector3& newPos) { this->points.insert(points.begin() + i, newPos); return *this; }
    BSpline& removePoint(int i) { this->points.erase(points.begin() + i); return *this; }

    void reset() { this->points.clear(); }

protected:
    float alpha = 1.f;  // alpha : 2 = very round, 1 = quite normal, 0.5 = almost linear
    std::vector<Vector3> points;
};

class CatmullRomSpline : public BSpline {

};

#endif // BSPLINE_H
