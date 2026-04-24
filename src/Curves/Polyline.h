#ifndef POLYLINE_H
#define POLYLINE_H

#include "Curves/Curve.h"

class Polyline : public Curve
{
public:
    Polyline();
    Polyline(const std::vector<Vector3>& points);

    virtual std::vector<Vector3> getPath(int numberOfPoints = -1) const override;
    virtual Vector3 getPoint(float x) const override;
    virtual Vector3 getPoint(float x, const Vector3& a, const Vector3& b) const override;
    virtual Vector3 getDerivative(float x, bool normalize = false) const override;
    virtual Vector3 getSecondDerivative(float x, bool normalize = false) const override;
    virtual float estimateClosestTime(const Vector3& pos, float epsilon = 1e-4, float nbChecksFactor = 2.f, float earlyExitThreshold = 1e-3) const override;
    virtual Vector3 estimateClosestPos(const Vector3& pos, bool useNativeShape = false, float epsilon = 1e-3) const override;
    virtual float estimateSqrDistanceFrom(const Vector3& pos, bool useNativeShape = false, float epsilon = 1e-3) const override;
    virtual float length() const override;

    size_t getIndex(int i) { return (i + numPoints()) % numPoints(); }
    virtual Polyline& setPoint(int i, const Vector3& newPos) override;

    virtual Polyline& resamplePoints(int newNbPoints = -1) override;

    virtual Polyline& reverseVertices() override;

    // virtual Polyline simplifyByRamerDouglasPeucker(float epsilon, Polyline subPolyline = Polyline()) override;

    virtual std::pair<Vector3, Vector3> AABBox() const override;

    Polyline& scale(float factor);
    Polyline& scale(const Vector3& factor);
    // Polyline scaled(float factor);
    // Polyline scaled(const Vector3& factor);

    virtual Polyline& translate(const Vector3& translation) override;

    virtual Polyline& removeDuplicates() override;

    size_t size() const { return numPoints(); }
    virtual size_t numPoints() const override;

    virtual Vector3& operator[](size_t i) override;
    virtual const Vector3& operator[](size_t i) const override;

    virtual std::string toString() const override;

    virtual Polyline& close() override;

protected:
    std::vector<Vector3> points;
};

#endif // POLYLINE_H
