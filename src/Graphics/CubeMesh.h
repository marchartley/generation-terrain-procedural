#ifndef CUBEMESH_H
#define CUBEMESH_H

#include "Graphics/Mesh.h"

class CubeMesh : public Mesh
{
public:
    CubeMesh();

    void update(std::vector<Vector3> positions = std::vector<Vector3>());
    void pushToBuffer(std::vector<Vector3> positions);
    void display(GLenum shape = 0);

    std::vector<Vector3> positions;
    std::vector<float> positionsFloat;

    static std::vector<Vector3> cubesVertices;
    static std::vector<Vector3> cubesEdgesVertices;

    static std::vector<Vector3> createTriangles(Vector3 minPos, Vector3 maxPos);
    static std::vector<Vector3> createTriangles(AABBox box);
};


#endif // CUBEMESH_H
