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
    virtual Vector3 getPoint(float x) const = 0; // alpha : 2 = very round, 1 = quite normal, 0.5 = almost linear
    virtual Vector3 getPoint(float x, const Vector3& a, const Vector3& b) const = 0;
    virtual Vector3 getDerivative(float x, bool normalize = false) const = 0;
    virtual Vector3 getSecondDerivative(float x, bool normalize = false) const = 0;
    virtual float estimateClosestTime(const Vector3& pos) const = 0;
    virtual Vector3 estimateClosestPos(const Vector3& pos) const = 0;
    virtual float estimateSqrDistanceFrom(const Vector3& pos) const = 0;
    virtual float estimateDistanceFrom(const Vector3& pos) const;
    virtual float estimateSignedDistanceFrom(const Vector3& pos) const;
    virtual float length() const = 0;
    virtual size_t timeToLowerIndex(float t) const;

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

    size_t size() const;
    virtual size_t numPoints() const = 0;
    virtual size_t numSegments() const;
    virtual size_t pointIndex(int index) const;
    bool empty() const;

    virtual Vector3& operator[](size_t i) = 0;
    virtual const Vector3& operator[](size_t i) const = 0;

    virtual std::string toString() const = 0;

    friend std::ostream& operator<<(std::ostream& io, const Curve& curve);
    friend std::ostream& operator<<(std::ostream& io, std::shared_ptr<Curve> curve);

    bool isClosed() const;
    Curve& setClosed(bool close);
    virtual Curve& close();

    virtual Curve& reset() = 0;


protected:
    bool closed = false;
};


#endif // CURVE_H
