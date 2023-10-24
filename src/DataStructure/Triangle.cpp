#include "Triangle.h"

//Triangle::Triangle()
//{

//}



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
