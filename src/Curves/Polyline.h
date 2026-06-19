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
    Vector3 getNormal(float x, Vector3 forcedUp = Vector3::invalid) const override;
    float estimateClosestTime(const Vector3& pos) const override;
    Vector3 estimateClosestPos(const Vector3& pos) const override;
    float estimateSqrDistanceFrom(const Vector3& pos) const override;
    float length() const override;

    size_t getIndex(int i);
    Polyline& setPoint(int i, const Vector3& newPos) override;

    std::vector<Vector3> getPoints() const;
    Vector3& get(int i) override;
    Vector3 get(int i) const override;

    Polyline& resamplePoints(int newNbPoints = -1) override;

    Polyline& reverseVertices() override;

    AABBox getAABBox() const override;

    using Curve::scale;
    Polyline& scale(const Vector3& factor) override;

    Polyline& translate(const Vector3& translation) override;

    Polyline& removeDuplicates() override;

    size_t size() const;
    size_t numPoints() const override;

    std::vector<Vector3>::const_iterator begin() const override;
    std::vector<Vector3>::const_iterator end() const override;
    std::vector<Vector3>::iterator begin() override;
    std::vector<Vector3>::iterator end() override;

    Vector3& operator[](size_t i) override;
    const Vector3& operator[](size_t i) const override;

    void addPoint(const Vector3& newPoint) override;
    Polyline& insertPoint(int i, const Vector3& newPos) override;
    Polyline& removePoint(int i) override;

    std::string toString() const override;

    Polyline& close() override;

    Polyline& reset();

    using Curve::slice;
    // std::vector<std::shared_ptr<Curve> > slice(const std::vector<float>& _ts) const;
    std::vector<std::shared_ptr<Curve>> slice(float t) const override;
    // std::vector<std::shared_ptr<Curve>> slice(std::vector<float> ts) const override;

protected:
    std::vector<Vector3> points;
};

#endif // POLYLINE_H
