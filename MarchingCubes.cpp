#include "MarchingCubes.h"

MarchingCubes::MarchingCubes()
{

}
MarchingCubes::MarchingCubes(VoxelChunk& grid)
    : MarchingCubes()
{
    map.clear();
    for (int x = 0; x < grid.sizeX; x++)
    {
        map.push_back(std::vector<std::vector<bool>>());
        for (int y = 0; y < grid.sizeY; y++)
        {
            map[x].push_back(std::vector<bool>());
            for (int z = 0; z < grid.height; z++)
            {
                map[x][y].push_back((bool)(grid.voxels[x][y][z]->type == TerrainTypes::AIR));
            }
        }
    }
}

void MarchingCubes::display()
{
    for (int x = 0; x < map.size() - 1; x++) {
        for (int y = 0; y < map[x].size() - 1; y++) {
            for (int z = 0; z < map[x][y].size() - 1; z++) {
//                bool vertices[8] = {map[x][y][z],
//                                    map[x+1][y][z],
//                                    map[x+1][y+1][z],
//                                    map[x][y+1][z],
//                                    map[x][y][z+1],
//                                    map[x+1][y][z+1],
//                                    map[x+1][y+1][z+1],
//                                    map[x][y+1][z+1]};
                Vertex vertices[8] = {Vertex(x, y, z, map[x][y][z] ? 1.0 : -1.0),
                                      Vertex(x+1, y, z, map[x+1][y][z] ? 1.0 : -1.0),
                                      Vertex(x+1, y+1, z, map[x+1][y+1][z] ? 1.0 : -1.0),
                                      Vertex(x, y+1, z, map[x][y+1][z] ? 1.0 : -1.0),
                                      Vertex(x, y, z+1, map[x][y][z+1] ? 1.0 : -1.0),
                                      Vertex(x+1, y, z+1, map[x+1][y][z+1] ? 1.0 : -1.0),
                                      Vertex(x+1, y+1, z+1, map[x+1][y+1][z+1] ? 1.0 : -1.0),
                                      Vertex(x, y+1, z+1, map[x][y+1][z+1] ? 1.0 : -1.0)
                                     };
                if (x == 0)
                {
                    vertices[0].isosurface = -1.0;
                    vertices[3].isosurface = -1.0;
                    vertices[4].isosurface = -1.0;
                    vertices[7].isosurface = -1.0;
                }if (x == map.size() - 1)
                {
                    vertices[1].isosurface = -1.0;
                    vertices[2].isosurface = -1.0;
                    vertices[5].isosurface = -1.0;
                    vertices[6].isosurface = -1.0;
                }if (y == 0)
                {
                    vertices[0].isosurface = -1.0;
                    vertices[1].isosurface = -1.0;
                    vertices[4].isosurface = -1.0;
                    vertices[5].isosurface = -1.0;
                }if (y == map[x].size() - 1)
                {
                    vertices[2].isosurface = -1.0;
                    vertices[3].isosurface = -1.0;
                    vertices[6].isosurface = -1.0;
                    vertices[7].isosurface = -1.0;
                }if (z == 0)
                {
                    vertices[0].isosurface = -1.0;
                    vertices[1].isosurface = -1.0;
                    vertices[2].isosurface = -1.0;
                    vertices[3].isosurface = -1.0;
                }if (z == map[x][y].size() - 1)
                {
                    vertices[4].isosurface = -1.0;
                    vertices[5].isosurface = -1.0;
                    vertices[6].isosurface = -1.0;
                    vertices[7].isosurface = -1.0;
                }
                int cube_index = 0;
                for (int i = 0; i < 8; i++){
                    if (vertices[i].isosurface > 0.0)
                        cube_index ^= 1 << i;
                }
                int* edgesForTriangles = MarchingCubes::triangleTable[cube_index];
                glBegin(GL_TRIANGLES);
                Vertex originalVertex;
                Vertex firstVertex;
                Vertex secondVertex;
                for (int i = 0; i < 16; i++) {
                    if (edgesForTriangles[i] == -1)
                        continue;
                    Vertex& v1 = vertices[MarchingCubes::edgeToCorner[edgesForTriangles[i]][0]];
                    Vertex& v2 = vertices[MarchingCubes::edgeToCorner[edgesForTriangles[i]][1]];
                    Vertex midpoint = (v1 + v2);
                    midpoint /= 2.0;
                    if (i % 3 == 0)
                        originalVertex = midpoint;
                    else if (i % 3 == 1)
                        firstVertex = midpoint;
                    else {
                        secondVertex = midpoint;
                        Vector3 normal = (firstVertex - originalVertex).cross((secondVertex - originalVertex)).normalize();
//                        glColor3f(abs(normal.x), abs(normal.y), abs(normal.z));
                        /*
                        std::cout << normal << " = " << firstVertex << " - " << originalVertex << std::endl;
                        std::cout << normal.x << ", " << 1-normal.x << ", 0.0" << std::endl;*/
                    }
                    glVertex3f(midpoint.x, midpoint.y, midpoint.z);
                }
                glEnd();
            }
        }
    }
}

void MarchingCubes::displayVoxel(Voxel &v) {
}
