#ifndef CURVE_H
#define CURVE_H

#include <vector>
#include "DataStructure/Vector3.h"

#define CLONE_FUNCTION(CLASS_NAME) \
virtual CLASS_NAME* clone() const { return new CLASS_NAME(*this); }

class Curve
{
public:
    Curve() = default;
    Curve(const Curve& copy) = default;
    virtual ~Curve() = default;

    virtual Curve* clone() const = 0;
    virtual std::vector<Vector3> getPath(int numberOfPoints = -1) const = 0;
    virtual Vector3 getPoint(float x) const = 0;
    Vector3 linear(float x, const Vector3& a, const Vector3& b) const;
    virtual Vector3 getDerivative(float x, bool normalize = false) const = 0;
    virtual Vector3 getSecondDerivative(float x, bool normalize = false) const = 0;
    virtual float estimateClosestTime(const Vector3& pos) const = 0;
    virtual Vector3 estimateClosestPos(const Vector3& pos) const = 0;
    virtual float estimateSqrDistanceFrom(const Vector3& pos) const = 0;
    virtual float estimateDistanceFrom(const Vector3& pos) const;
    virtual float estimateSignedDistanceFrom(const Vector3& pos) const;
    virtual float length() const = 0;
    virtual size_t timeToLowerIndex(float t) const;
    virtual Vector3& get(int i) = 0; // Get point by index instead of time
    virtual Vector3 get(int i) const { return get(i); }

    virtual Curve& setPoint(int i, const Vector3& newPos) = 0;

    virtual Curve& resamplePoints(int newNbPoints = -1) = 0;

    virtual Curve& reverseVertices() = 0;

    virtual Vector3 getDirection(float x) const;
    virtual Vector3 getNormal(float x) const;
    virtual Vector3 getBinormal(float x) const;
    virtual float getCurvature(float x) const;

    // virtual Curve simplifyByRamerDouglasPeucker(float epsilon, Curve subcurve = Curve()) = 0;

    virtual std::pair<Vector3, Vector3> AABBox() const = 0;
    virtual Vector3 containingBoxSize() const;

    virtual Curve& scale(float factor);
    virtual Curve& scale(const Vector3& factor) = 0;
    // Curve scaled(float factor);
    // Curve scaled(const Vector3& factor);

    virtual Curve& translate(const Vector3& translation) = 0;

    virtual Curve& removeDuplicates() = 0;

    virtual Vector3 center() const;

    size_t size() const;
    virtual size_t numPoints() const = 0;
    virtual size_t numSegments() const;
    virtual size_t pointIndex(int index) const;
    virtual bool empty() const;

    virtual std::vector<Vector3>::const_iterator begin() const = 0;
    virtual std::vector<Vector3>::const_iterator end() const = 0;
    virtual std::vector<Vector3>::iterator begin() = 0;
    virtual std::vector<Vector3>::iterator end() = 0;

    inline Vector3 front() const { return (empty() ? Vector3::invalid : get(0)); }
    inline Vector3 back() const { return (empty() ? Vector3::invalid : get(-1)); }

    virtual Vector3& operator[](size_t i) = 0;
    virtual const Vector3& operator[](size_t i) const = 0;

    virtual void addPoint(const Vector3& newPoint) = 0;
    virtual Curve& insertPoint(int i, const Vector3& newPos) = 0;
    virtual Curve& removePoint(int i) = 0;

    virtual std::string toString() const = 0;

    friend std::ostream& operator<<(std::ostream& io, const Curve& curve);
    friend std::ostream& operator<<(std::ostream& io, std::shared_ptr<Curve> curve);

    bool isClosed() const;
    Curve& setClosed(bool close);
    virtual Curve& close();

    virtual Curve& reset() = 0;

    virtual std::vector<std::pair<size_t, size_t>> checkAutointersections() const;


protected:
    bool closed = false;
};


#endif // CURVE_H
