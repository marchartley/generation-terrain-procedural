#ifndef POLYLINE_H
#define POLYLINE_H

#include "Curves/Curve.h"

class Polyline : public Curve
{
public:
    Polyline();
    Polyline(const std::vector<Vector3>& points);

    CLONE_FUNCTION(Polyline)
    std::vector<Vector3> getPath(int numberOfPoints = -1) const override;
    Vector3 getPoint(float x) const override;
    Vector3 getPoint(float x, const Vector3& a, const Vector3& b) const override;
    Vector3 getDerivative(float x, bool normalize = false) const override;
    Vector3 getSecondDerivative(float x, bool normalize = false) const override;
    float estimateClosestTime(const Vector3& pos) const override;
    Vector3 estimateClosestPos(const Vector3& pos) const override;
    float estimateSqrDistanceFrom(const Vector3& pos) const override;
    float length() const override;

    size_t getIndex(int i) { return (i + numPoints()) % numPoints(); }
    Polyline& setPoint(int i, const Vector3& newPos) override;

    Polyline& resamplePoints(int newNbPoints = -1) override;

    Polyline& reverseVertices() override;

    // Polyline simplifyByRamerDouglasPeucker(float epsilon, Polyline subPolyline = Polyline()) override;

    std::pair<Vector3, Vector3> AABBox() const override;

    using Curve::scale;
    Polyline& scale(const Vector3& factor) override;
    // Polyline scaled(float factor);
    // Polyline scaled(const Vector3& factor);

    Polyline& translate(const Vector3& translation) override;

    Polyline& removeDuplicates() override;

    size_t size() const { return numPoints(); }
    size_t numPoints() const override;

    Vector3& operator[](size_t i) override;
    const Vector3& operator[](size_t i) const override;

    std::string toString() const override;

    Polyline& close() override;

    Polyline& reset() { points.clear(); return *this; }

    std::vector<Vector3> getPoints() const { return this->points; }

protected:
    std::vector<Vector3> points;
};

#endif // POLYLINE_H
