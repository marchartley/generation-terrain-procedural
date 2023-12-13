#include "Triangle.h"

//Triangle::Triangle()
//{

//}



Triangle::Triangle()
    : Triangle(Vector3(), Vector3(), Vector3())
{}

Triangle::Triangle(const Vector3 &p1, const Vector3 &p2, const Vector3 &p3) {
    vertices[0] = p1;
    vertices[1] = p2;
    vertices[2] = p3;

    // Compute normal
    Vector3 edge1 = p2 - p1;
    Vector3 edge2 = p3 - p1;
    normal = edge1.cross(edge2).normalized();

    d = p1.dot(normal);
}

Triangle::Triangle(std::vector<Vector3> tri)
    : Triangle(tri[0], tri[1], tri[2])
{}

const Vector3 &Triangle::operator[](size_t i) const { return this->vertices[i]; }

Vector3 &Triangle::operator[](size_t i) { return this->vertices[i]; }

Triangle& Triangle::enlarge(float addedSize)
{
    Vector3 center = (this->vertices[0] + this->vertices[1] + this->vertices[2]) / 3.f;
    for (auto& v : vertices) {
        v += (v - center).setMag(addedSize);
    }
    return *this;
}

std::vector<Triangle> Triangle::vectorsToTriangles(const std::vector<std::vector<Vector3> > &vectors)
{
    std::vector<Triangle> tris(vectors.size());
    for (size_t i = 0; i < vectors.size(); i++)
        tris[i] = Triangle(vectors[i]);
    return tris;
}

std::vector<std::vector<Vector3> > Triangle::trianglesToVectors(const std::vector<Triangle> &triangles)
{
    std::vector<std::vector<Vector3>> tris(triangles.size());
    for (size_t i = 0; i < triangles.size(); i++)
        tris[i] = {triangles[i][0], triangles[i][1], triangles[i][2]};
    return tris;
}
