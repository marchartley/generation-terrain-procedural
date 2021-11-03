#include "MarchingCubes.h"



MarchingCubes::MarchingCubes()
{

}
MarchingCubes::MarchingCubes(VoxelChunk& grid)
    : MarchingCubes()
{
    this->grid = &grid;
    map.clear();
    for (int x = 0; x < grid.sizeX; x++)
    {
        map.push_back(std::vector<std::vector<bool>>());
        for (int y = 0; y < grid.sizeY; y++)
        {
            map[x].push_back(std::vector<bool>());
            for (int z = 0; z < grid.height; z++)
            {
                map[x][y].push_back((bool)(grid.voxels[x][y][z]->type == TerrainTypes::DIRT));
            }
        }
    }

    bool addedLeft = false, addedRight = false, addedFront = false, addedBack = false;
    if (this->grid->neighboring_chunks.find(LEFT) != this->grid->neighboring_chunks.end()) {
        VoxelChunk* n = this->grid->neighboring_chunks[LEFT];
        map.insert(map.begin(), std::vector<std::vector<bool>>());
        for (int y = 0; y < this->grid->sizeY; y++) {
            map[0].push_back(std::vector<bool>());
            for (int z = 0; z < grid.height; z++)
            {
                map[0][y].push_back((bool)(n->voxels[n->sizeX - 1][y][z]->type == TerrainTypes::DIRT));
            }
        }
        addedLeft = true;
    }
//    if (this->grid->neighboring_chunks.find(RIGHT) != this->grid->neighboring_chunks.end()) {
//        VoxelChunk* n = this->grid->neighboring_chunks[RIGHT];
//        map.push_back(std::vector<std::vector<bool>>());
//        int index = map.size() - 1;
//        for (int y = 0; y < this->grid->sizeY; y++) {
//            map[index].push_back(std::vector<bool>());
//            for (int z = 0; z < grid.height; z++)
//            {
//                map[index][y].push_back((bool)(n->voxels[0][y][z]->type == TerrainTypes::DIRT));
//            }
//        }
//        addedRight = true;
//    }
//    if (this->grid->neighboring_chunks.find(BACK) != this->grid->neighboring_chunks.end()) {
//        VoxelChunk* n = this->grid->neighboring_chunks[BACK];
//        int offset = addedLeft ? 1 : 0;
//        for (int x = 0; x < this->grid->sizeX; x++) {
//            map[x + offset].push_back(std::vector<bool>());
//            int y = this->map[offset].size() - 1;
//            for (int z = 0; z < grid.height; z++)
//            {
//                map[x + offset][y].push_back((bool)(n->voxels[x][0][z]->type == TerrainTypes::DIRT));
//            }
//        }
//        addedBack = true;
//    }
    if (this->grid->neighboring_chunks.find(FRONT) != this->grid->neighboring_chunks.end()) {
        VoxelChunk* n = this->grid->neighboring_chunks[FRONT];
        int offset = addedLeft ? 1 : 0;
        for (int x = 0; x < this->grid->sizeX; x++) {
            map[x + offset].insert(map[x + offset].begin(), std::vector<bool>());
            for (int z = 0; z < grid.height; z++)
            {
                map[x + offset][0].push_back((bool)(n->voxels[x][n->voxels[x].size() - 1][z]->type == TerrainTypes::DIRT));
            }
        }
        addedFront = true;
    }

    if (addedLeft) {
        if (addedBack) {
            VoxelChunk* n = this->grid->neighboring_chunks[LEFT]->neighboring_chunks[BACK];
            map[0].push_back(std::vector<bool>());
            for (int z = 0; z < grid.height; z++)
                map[0][map[0].size() - 1].push_back((bool)(n->voxels[n->voxels.size() - 1][0][z]->type == TerrainTypes::DIRT));
        }
        if (addedFront) {
            VoxelChunk* n = this->grid->neighboring_chunks[LEFT]->neighboring_chunks[FRONT];
            map[0].insert(map[0].begin(), std::vector<bool>());
            for (int z = 0; z < grid.height; z++)
                map[0][0].push_back((bool)(n->voxels[n->voxels.size() - 1][n->voxels[0].size() - 1][z]->type == TerrainTypes::DIRT));
        }
    }
//    if (addedRight) {
//        if (addedBack) {
//            VoxelChunk* n = this->grid->neighboring_chunks[RIGHT]->neighboring_chunks[BACK];
//            map[map.size() - 1].insert(map[map.size() - 1].begin(), std::vector<bool>());
//            for (int z = 0; z < grid.height; z++)
//                map[map.size() - 1][0].push_back((bool)(n->voxels[0][0][z]->type == TerrainTypes::DIRT));
//        }
//        if (addedFront) {
//            VoxelChunk* n = this->grid->neighboring_chunks[RIGHT]->neighboring_chunks[FRONT];
//            map[map.size() - 1].push_back(std::vector<bool>());
//            for (int z = 0; z < grid.height; z++)
//                map[map.size() - 1][map[0].size() - 1].push_back((bool)(n->voxels[0][n->voxels[0].size() - 1][z]->type == TerrainTypes::DIRT));
//        }
//    }
//    for (int h = 0; h < this->map[0][0].size(); h++) {
//        for (int y = 0; y < this->map[0].size(); y++) {
//            for (int x = 0; x < this->map.size(); x++) {
//                std::cout << (this->map[x][y][h] ? "T" : "F") << "|";
//            }
//            std::cout << std::endl;
//        }
//        std::cout << std::endl;
//    }
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
                if (x == 0 && this->grid->x == 0)
                {
                    vertices[0].isosurface = -1.0;
                    vertices[3].isosurface = -1.0;
                    vertices[4].isosurface = -1.0;
                    vertices[7].isosurface = -1.0;
                } else if (x == map.size() - (this->grid->x > 0 ? 2 : 1) && this->grid->lastChunkOnX)
                {
                    vertices[1].isosurface = -1.0;
                    vertices[2].isosurface = -1.0;
                    vertices[5].isosurface = -1.0;
                    vertices[6].isosurface = -1.0;
                }
                if (y == 0 && this->grid->y == 0)
                {
                    vertices[0].isosurface = -1.0;
                    vertices[1].isosurface = -1.0;
                    vertices[4].isosurface = -1.0;
                    vertices[5].isosurface = -1.0;
                }if (y == map[x].size() - (this->grid->y > 0 ? 2 : 1) && this->grid->lastChunkOnY)
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
                }if (z == map[x][y].size() - 2)
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
                    Vertex midpoint = (v1 + v2) / 2.0;
//                    std::cout << "V1 " << v1 << " V2 " << v2 << std::endl << "V2 - V1 " << (v2-v1) << " / " << (v2.isosurface - v1.isosurface) << " = " << (v2 - v1) / (v1.isosurface - v2.isosurface) << std::endl;
//                    std::cout << " - V1 = " << v1 - ((v2 - v1) / (v2.isosurface - v1.isosurface)) << std::endl;
//                    Vertex midpoint = v1 - ((v2 - v1) / (v2.isosurface - v1.isosurface));
//                    midpoint /= 2.0;
//                    Vertex midpoint = v1;
//                    midpoint += (v1 - v2)/2.0f;
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
