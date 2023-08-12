#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "DataStructure/Vector3.h"


class Triangle {
public:
    std::vector<Vector3> vertices = std::vector<Vector3>(3);
    Vector3 normal;
    float d;

    Triangle()
        : Triangle(Vector3(), Vector3(), Vector3())
    {}
    Triangle(const Vector3& p1, const Vector3& p2, const Vector3& p3) {
        vertices[0] = p1;
        vertices[1] = p2;
        vertices[2] = p3;

        // Compute normal
        Vector3 edge1 = p2 - p1;
        Vector3 edge2 = p3 - p1;
        normal = edge1.cross(edge2).normalized();

        d = p1.dot(normal);
    }
    Triangle(std::vector<Vector3> tri)
        : Triangle(tri[0], tri[1], tri[2])
    {}

    const Vector3& operator[](size_t i) const { return this->vertices[i]; }
    Vector3& operator[](size_t i) { return this->vertices[i]; }

    static std::vector<Triangle> vectorsToTriangles(const std::vector<std::vector<Vector3>>& vectors);
    static std::vector<std::vector<Vector3>> trianglesToVectors(const std::vector<Triangle>& triangles);
};



#endif // TRIANGLE_H
