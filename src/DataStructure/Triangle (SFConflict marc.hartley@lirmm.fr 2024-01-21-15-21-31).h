#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "DataStructure/Vector3.h"


class Triangle {
public:
    std::vector<Vector3> vertices = std::vector<Vector3>(3);
    Vector3 normal;
    float d;

    Triangle();
    Triangle(const Vector3& p1, const Vector3& p2, const Vector3& p3);
    Triangle(std::vector<Vector3> tri);

    const Vector3& operator[](size_t i) const;
    Vector3& operator[](size_t i);

    Triangle &enlarge(float addedSize);

    static std::vector<Triangle> vectorsToTriangles(const std::vector<std::vector<Vector3>>& vectors);
    static std::vector<std::vector<Vector3>> trianglesToVectors(const std::vector<Triangle>& triangles);
};



#endif // TRIANGLE_H
