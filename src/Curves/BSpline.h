#ifndef BSPLINE_H
#define BSPLINE_H

#include "Curves/Curve.h"

class BSpline : public Curve
{
public:
    BSpline();
    BSpline(const std::vector<Vector3>& points, bool clamped = true);
    BSpline(const std::vector<Vector3>& points, const std::vector<float>& knots);
    BSpline(const Curve& curve);

    CLONE_FUNCTION(BSpline)
    std::vector<Vector3> getPath(int numberOfPoints = -1) const override;
    Vector3 getPoint(float x) const override;
    Vector3 getDerivative(float x, bool normalize = false) const override;
    Vector3 getSecondDerivative(float x, bool normalize = false) const override;
    float estimateClosestTime(const Vector3& pos) const override;
    Vector3 estimateClosestPos(const Vector3& pos) const override;
    float estimateSqrDistanceFrom(const Vector3& pos) const override;
    float length() const override;

    size_t getIndex(int i);
    BSpline& setPoint(int i, const Vector3& newPos) override;

    std::vector<Vector3> getPoints() const;
    Vector3& get(int i) override;
    Vector3 get(int i) const override;

    size_t numSegments() const override;

    BSpline& resamplePoints(int newNbPoints = -1) override;

    BSpline& reverseVertices() override;

    std::pair<Vector3, Vector3> AABBox() const override;

    using Curve::scale;
    BSpline& scale(const Vector3& factor) override;

    BSpline& translate(const Vector3& translation) override;

    BSpline& removeDuplicates() override;

    size_t size() const;
    size_t numPoints() const override;

    std::vector<Vector3>::const_iterator begin() const override;
    std::vector<Vector3>::const_iterator end() const override;
    std::vector<Vector3>::iterator begin() override;
    std::vector<Vector3>::iterator end() override;

    Vector3& operator[](size_t i) override;
    const Vector3& operator[](size_t i) const override;

    void addPoint(const Vector3& newPoint) override;
    BSpline& insertPoint(int i, const Vector3& newPos) override;
    BSpline& removePoint(int i) override;

    std::string toString() const override;

    BSpline& close() override;

    BSpline& reset();

    std::vector<std::shared_ptr<Curve>> slice(const std::vector<float>& _ts) const override;
    std::vector<std::shared_ptr<Curve>> slice(float t) const override;

    BSpline& setDegree(int newDegree);
    int getDegree() const { return degree; }

    std::vector<float> getKnots() const { return knots; }

    std::vector<float> uniqueInternalKnots() const;
    std::vector<float> uniqueInternalKnots(float minKnot, float maxKnot) const;

    int knotMultiplicity(float u) const;

    size_t insertKnot(float u);
    BSpline& removeKnotAtKnotIndex(int knotIndex, float tolerance = 1e-4);
    BSpline& removeKnot(int pointIndex, float tolerance = 1e-4);

    BSpline generateFakeClosedCurve() const;
protected:
    std::vector<Vector3> points;
    std::vector<float> knots;

    int degree = 3;

    std::vector<Vector3> solveLinearSystem(std::vector<std::vector<float>> A, std::vector<Vector3> b);
    float basis(int i, int p, float u, const std::vector<float>& U);
    std::vector<float> makeUniformClampedKnots(int numCtrlPoints, int degree);

    static float u_from_t(float t, const std::vector<float>& U, int nbPoints, int degree);
    float u_from_t(float t) const;

};

#endif // BSPLINE_H
