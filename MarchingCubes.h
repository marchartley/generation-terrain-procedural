#ifndef MARCHINGCUBES_H
#define MARCHINGCUBES_H

#include "Voxel.h"
#include "Mesh.h"

class MarchingCubes
{
public:
    static int cubeEdges[256];
    static int triangleTable[256][16];
    static int edgeToCorner[12][2];
};

#endif // MARCHINGCUBES_H
