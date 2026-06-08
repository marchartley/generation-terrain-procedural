#ifndef BEZIERCURVE_H
#define BEZIERCURVE_H

#include "Curves/Curve.h"

class BezierCurve : public Curve
{
public:
    BezierCurve();
    BezierCurve(const std::vector<Vector3>& points);
    BezierCurve(const std::vector<Vector3>& points, const std::vector<Vector3>& handles);
    // BezierCurve(const BezierCurve& curve);
    ~BezierCurve() = default;

    CLONE_FUNCTION(BezierCurve)
    std::vector<Vector3> getPath(int numberOfPoints = -1) const override;
    Vector3 getPoint(float x) const override;
    Vector3 getDerivative(float x, bool normalize = false) const override;
    Vector3 getSecondDerivative(float x, bool normalize = false) const override;
    float estimateClosestTime(const Vector3& pos) const override;
    Vector3 estimateClosestPos(const Vector3& pos) const override;
    float estimateSqrDistanceFrom(const Vector3& pos) const override;
    float length() const override;

    Vector3& get(int i) override;
    Vector3 get(int i) const override;

    size_t getIndex(int i);
    BezierCurve& setPoint(int i, const Vector3& newPos) override;

    BezierCurve& resamplePoints(int newNbPoints = -1) override;

    BezierCurve& reverseVertices() override;

    // BezierCurve simplifyByRamerDouglasPeucker(float epsilon, BezierCurve subBezierCurve = BezierCurve()) override;

    std::pair<Vector3, Vector3> AABBox() const override;

    using Curve::scale;
    BezierCurve& scale(const Vector3& factor) override;
    // BezierCurve scaled(float factor);
    // BezierCurve scaled(const Vector3& factor);

    BezierCurve& translate(const Vector3& translation) override;

    BezierCurve& removeDuplicates() override;

    size_t size() const;
    size_t numPoints() const override;

    std::vector<Vector3>::const_iterator begin() const override;
    std::vector<Vector3>::const_iterator end() const override;
    std::vector<Vector3>::iterator begin() override;
    std::vector<Vector3>::iterator end() override;

    Vector3& operator[](size_t i) override;
    const Vector3& operator[](size_t i) const override;

    void addPoint(const Vector3& newPoint) override;
    BezierCurve& insertPoint(int i, const Vector3& newPos) override;
    BezierCurve& removePoint(int i) override;

    std::string toString() const override;

    BezierCurve& close() override;

    BezierCurve& reset();

    using Curve::pointIndex;
    inline static size_t pointIndex(int index, int nbPoints, bool closed);
    inline static size_t handleIndex(int index, int nbPoints, bool closed);
    inline static size_t handleIn(int pointIdx, int nbPoints, bool closed);
    inline static size_t handleOut(int pointIdx, int nbPoints, bool closed);

    inline size_t handleIndex(int index) const;
    inline size_t handleIn(int pointIdx) const;
    inline size_t handleOut(int pointIdx) const;

    inline size_t pointIndexFromHandleIndex(int handleIdx) const;
    inline Vector3 handlePos(int handleIdx) const;

    std::vector<Vector3> getPoints() const;
    std::vector<Vector3> getHandles() const;

    BezierCurve& autosmooth(int pointIndex);
    BezierCurve& autosmooth();

    std::vector<std::shared_ptr<Curve>> slice(float t) const override;
    std::vector<std::shared_ptr<Curve>> slice(const std::vector<float>& ts) const override;


    static Vector3 linearBezier(const Vector3& P0, const Vector3& P1, float t);
    static Vector3 quadraticBezier(const Vector3& P0, const Vector3& P1, const Vector3& P2, float t);
    static Vector3 cubicBezier(const Vector3& P0, const Vector3& P1, const Vector3& P2, const Vector3& P3, float t);
    static Vector3 cubicBezierDerivative(const Vector3& P0, const Vector3& P1, const Vector3& P2, const Vector3& P3, float t);
    static Vector3 cubicBezierSecondDerivative(const Vector3& P0, const Vector3& P1, const Vector3& P2, const Vector3& P3, float t);


// protected:
    std::vector<Vector3> points;
    std::vector<Vector3> handles;
};

#endif // BEZIERCURVE_H
