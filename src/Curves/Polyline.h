#ifndef POLYLINE_H
#define POLYLINE_H

#include "Curves/Curve.h"

class Polyline : public Curve
{
public:
    Polyline();
    Polyline(const std::vector<Vector3>& points);
    Polyline(const Curve& curve);

    CLONE_FUNCTION(Polyline)
    std::vector<Vector3> getPath(int numberOfPoints = -1) const override;
    Vector3 getPoint(float x) const override;
    Vector3 getDerivative(float x, bool normalize = false) const override;
    Vector3 getSecondDerivative(float x, bool normalize = false) const override;
    float estimateClosestTime(const Vector3& pos) const override;
    Vector3 estimateClosestPos(const Vector3& pos) const override;
    float estimateSqrDistanceFrom(const Vector3& pos) const override;
    float length() const override;

    size_t getIndex(int i) { return (i + numPoints()) % numPoints(); }
    Polyline& setPoint(int i, const Vector3& newPos) override;

    std::vector<Vector3> getPoints() const { return this->points; }
    Vector3& get(int i) override { return this->points[pointIndex(i)]; }
    Vector3 get(int i) const override { return this->points[pointIndex(i)]; }

    Polyline& resamplePoints(int newNbPoints = -1) override;

    Polyline& reverseVertices() override;

    std::pair<Vector3, Vector3> AABBox() const override;

    using Curve::scale;
    Polyline& scale(const Vector3& factor) override;

    Polyline& translate(const Vector3& translation) override;

    Polyline& removeDuplicates() override;

    inline size_t size() const { return numPoints(); }
    size_t numPoints() const override;

    std::vector<Vector3>::const_iterator begin() const override { return points.begin(); }
    std::vector<Vector3>::const_iterator end() const override { return points.end(); }
    std::vector<Vector3>::iterator begin() override { return points.begin(); }
    std::vector<Vector3>::iterator end() override { return points.end(); }

    Vector3& operator[](size_t i) override;
    const Vector3& operator[](size_t i) const override;

    void addPoint(const Vector3& newPoint) override { this->points.push_back(newPoint); }
    Polyline& insertPoint(int i, const Vector3& newPos) override { this->points.insert(points.begin() + i, newPos); return *this; }
    Polyline& removePoint(int i) override { this->points.erase(points.begin() + i); return *this; }

    std::string toString() const override;

    Polyline& close() override;

    Polyline& reset() { points.clear(); return *this; }

protected:
    std::vector<Vector3> points;
};

#endif // POLYLINE_H
