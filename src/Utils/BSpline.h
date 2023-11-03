#ifndef BSPLINE_H
#define BSPLINE_H

#include <vector>
#include "DataStructure/Vector3.h"

class BSpline
{
public:
    BSpline();
    BSpline(int numberOfPoints);
    BSpline(std::vector<Vector3> points);
    BSpline(std::vector<BSpline> subsplines);

    std::vector<Vector3> getPath(int numberOfPoints) const;
    Vector3 getPoint(float x) const;
    Vector3 getPoint(float x, const Vector3& a, const Vector3& b) const;
    Vector3 getDerivative(float x) const;
    Vector3 getSecondDerivative(float x) const;
    float estimateClosestTime(const Vector3& pos, float epsilon = 1e-5, float nbChecksFactor = 4.f, float earlyExitThreshold = 1e-3) const;
    Vector3 estimateClosestPos(const Vector3& pos, float epsilon = 1e-3) const;
    float estimateSqrDistanceFrom(const Vector3& pos, float epsilon = 1e-3) const;
    float estimateDistanceFrom(const Vector3& pos, float epsilon = 1e-3) const;
    float estimateSignedDistanceFrom(const Vector3& pos, float epsilon = 1e-3) const;
    float length() const;

    BSpline& resamplePoints(int newNbPoints = -1);

    size_t nextID(int i) { return (i + 1 + this->points.size()) % this->points.size(); }

    operator bool() const { return (this->points.size() > 0); };

    std::tuple<Vector3, Vector3, Vector3> getFrenetFrame(float x) const;
    Vector3 getFrenetDirection(float x) const;
    Vector3 getFrenetNormal(float x) const;
    Vector3 getFrenetBinormal(float x) const;

    Vector3 getCenterCircle(float x) const;
    Vector3 getDirection(float x) const;
    Vector3 getNormal(float x) const;
    Vector3 getBinormal(float x) const;
    float getCurvature(float x) const;

    Vector3 center() const;

    BSpline& close();

    Vector3 getCatmullPoint(float x) const;

    BSpline simplifyByRamerDouglasPeucker(float epsilon, BSpline subspline = BSpline());

    std::tuple<Vector3, Vector3> AABBox() const;
    Vector3 containingBoxSize() const;

//    BSpline& grow(float increase);
//    BSpline& shrink(float decrease);

    BSpline computeConvexHull() const;

    BSpline& translate(const Vector3& translation);

    virtual BSpline& removeDuplicates();

    std::string toString() const;

    friend std::ostream& operator<<(std::ostream& io, const BSpline& s);
    friend std::ostream& operator<<(std::ostream& io, std::shared_ptr<BSpline> s);

    std::vector<Vector3> points;
    bool closed = false;


    auto begin() const { return points.begin(); }
    auto end() const { return points.end(); }
    auto begin() { return points.begin(); }
    auto end() { return points.end(); }
    std::size_t size() const { return end() - begin(); }
    bool empty() const { return begin() == end(); }

    Vector3& operator[](size_t i);
    const Vector3& operator[](size_t i) const;

    std::string display1DPlot(int sizeX, int sizeY) const;
};


#include "Utils/json.h"
nlohmann::json bspline_to_json(const BSpline& spline);
BSpline json_to_bspline(nlohmann::json json);

#endif // BSPLINE_H
