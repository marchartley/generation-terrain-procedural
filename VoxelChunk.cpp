#include "VoxelChunk.h"

VoxelChunk::VoxelChunk(int _x, int _y, int _sizeX, int _sizeY, int _height, std::vector<std::vector<std::vector< TerrainTypes > > > _data, VoxelGrid* parent)
    : data(_data), height(_height), x(_x), y(_y), sizeX(_sizeX), sizeY(_sizeY), parent(parent) {
    iso_data.clear();
    for (int x = 0; x < sizeX; x++) {
        iso_data.push_back(std::vector<std::vector<float> >());
        for (int y = 0; y < sizeY; y++) {
            iso_data[x].push_back(std::vector<float>());
            for(int z = 0; z < height; z++) {
                iso_data[x][y].push_back(data[x][y][z] == TerrainTypes::AIR ? -1.0 : 1.0);
            }
        }
    }
}
VoxelChunk::VoxelChunk(int x, int y, int sizeX, int sizeY, int height, std::vector<std::vector<std::vector< float > > > iso_data, VoxelGrid* parent)
    : iso_data(iso_data), height(height), x(x), y(y), sizeX(sizeX), sizeY(sizeY), parent(parent) {
    data.clear();
    for (int x = 0; x < sizeX; x++) {
        data.push_back(std::vector<std::vector<TerrainTypes> >());
        for (int y = 0; y < sizeY; y++) {
            data[x].push_back(std::vector<TerrainTypes>());
            for(int z = 0; z < height; z++) {
                data[x][y].push_back(iso_data[x][y][z] > 0.0 ? TerrainTypes::DIRT : TerrainTypes::AIR);
            }
        }
    }
}

VoxelChunk::VoxelChunk() : VoxelChunk(0, 0, 0, 0, 0, std::vector<std::vector<std::vector<TerrainTypes>>>(), nullptr)
{

}

void VoxelChunk::createMesh() {
    if (!needRemeshing)
        return;
    if (this->voxels.size() == 0) {
        this->voxels = std::vector<std::vector<std::vector<Voxel*>>>();
        for(int v_x = 0; v_x < sizeX; v_x++) {
            this->voxels.push_back(std::vector<std::vector<Voxel*>>());
            for(int v_y = 0; v_y < sizeY; v_y++) {
                this->voxels[v_x].push_back(std::vector<Voxel*>());
                for(int h = 0; h < height; h++) {
                    Voxel* v = new Voxel(v_x, v_y, h, this->data[v_x][v_y][h], 1.0, iso_data[v_x][v_y][h]);
                    //                    Voxel* v = new Voxel(v_x, v_y, h, this->data[v_x][v_y][h], 1.0);
                    v->parent = this;
                    this->voxels[v_x][v_y].push_back(v);
                }
            }
        }
    }

    for(int v_x = 0; v_x < sizeX; v_x++) {
        for(int v_y = 0; v_y < sizeY; v_y++) {
            for(int h = 0; h < height; h++) {
                Voxel* v = this->voxels[v_x][v_y][h];
                v->has_neighbors[LEFT] = v_x > 0 && this->data[v_x-1][v_y][h];
                v->has_neighbors[RIGHT] = v_x < sizeX-1 && this->data[v_x+1][v_y][h];
                v->has_neighbors[FRONT] = v_y > 0 && this->data[v_x][v_y-1][h];
                v->has_neighbors[BACK] = v_y < sizeY-1 && this->data[v_x][v_y+1][h];
                v->has_neighbors[BOTTOM] = h > 0 && this->data[v_x][v_y][h-1];
                v->has_neighbors[TOP] = h < height - 1 && this->data[v_x][v_y][h+1];

                if(v_x == 0 && this->neighboring_chunks.find(LEFT) != this->neighboring_chunks.end())
                {
                    VoxelChunk* n = this->neighboring_chunks[LEFT];
                    v->has_neighbors[LEFT] = (bool)(n->data[n->sizeX-1][v_y][h]);
                }
                if(v_y == 0 && this->neighboring_chunks.find(FRONT) != this->neighboring_chunks.end())
                {
                    VoxelChunk* n = this->neighboring_chunks[FRONT];
                    v->has_neighbors[FRONT] = (bool)(n->data[v_x][n->sizeY - 1][h]);
                }
                if(v_x == sizeX - 1 && this->neighboring_chunks.find(RIGHT) != this->neighboring_chunks.end())
                {
                    VoxelChunk* n = this->neighboring_chunks[RIGHT];
                    v->has_neighbors[RIGHT] = (bool)(n->data[0][v_y][h]);
                }
                if(v_y == sizeY - 1 && this->neighboring_chunks.find(BACK) != this->neighboring_chunks.end())
                {
                    VoxelChunk* n = this->neighboring_chunks[BACK];
                    v->has_neighbors[BACK] = (bool)(n->data[v_x][0][h]);
                }
            }
        }
    }
    needRemeshing = false;

//    MarchingCubes mc = MarchingCubes(*this);
//    mc.createMesh();
    this->mesh.fromArray(this->applyMarchingCubes()); // = mc.mesh;
    this->mesh.update();

}

std::vector<Vector3> VoxelChunk::applyMarchingCubes()
{
    std::vector<std::vector<std::vector<float> > > map;
    for (int x = 0; x < this->sizeX; x++)
    {
        map.push_back(std::vector<std::vector<float>>());
        for (int y = 0; y < this->sizeY; y++)
        {
            map[x].push_back(std::vector<float>());
            for (int z = 0; z < this->height; z++)
            {
                map[x][y].push_back(this->voxels[x][y][z]->getIsosurface());
            }
        }
    }

    bool addedLeft = false, addedRight = false, addedFront = false, addedBack = false;
    if (this->neighboring_chunks.find(LEFT) != this->neighboring_chunks.end()) {
        VoxelChunk* n = this->neighboring_chunks[LEFT];
        map.insert(map.begin(), std::vector<std::vector<float>>());
        for (int y = 0; y < this->sizeY; y++) {
            map[0].push_back(std::vector<float>());
            for (int z = 0; z < this->height; z++)
            {
                map[0][y].push_back(n->voxels[n->sizeX - 1][y][z]->getIsosurface());
            }
        }
        addedLeft = true;
    }
    if (this->neighboring_chunks.find(FRONT) != this->neighboring_chunks.end()) {
        VoxelChunk* n = this->neighboring_chunks[FRONT];
        int offset = addedLeft ? 1 : 0;
        for (int x = 0; x < this->sizeX; x++) {
            map[x + offset].insert(map[x + offset].begin(), std::vector<float>());
            for (int z = 0; z < this->height; z++)
            {
                map[x + offset][0].push_back(n->voxels[x][n->voxels[x].size() - 1][z]->getIsosurface());
            }
        }
        addedFront = true;
    }

    if (addedLeft) {
        if (addedBack) {
            VoxelChunk* n = this->neighboring_chunks[LEFT]->neighboring_chunks[BACK];
            map[0].push_back(std::vector<float>());
            for (int z = 0; z < this->height; z++)
                map[0][map[0].size() - 1].push_back(n->voxels[n->voxels.size() - 1][0][z]->getIsosurface());
        }
        if (addedFront) {
            VoxelChunk* n = this->neighboring_chunks[LEFT]->neighboring_chunks[FRONT];
            map[0].insert(map[0].begin(), std::vector<float>());
            for (int z = 0; z < this->height; z++)
                map[0][0].push_back(n->voxels[n->voxels.size() - 1][n->voxels[0].size() - 1][z]->getIsosurface());
        }
    }

    float isolevel = 0.0;
    std::vector<Vector3> vertexArray;
    for (unsigned int x = 0; x < map.size() - 1; x++) {
        for (unsigned int y = 0; y < map[x].size() - 1; y++) {
            for (unsigned int z = 0; z < map[x][y].size() - 1; z++) {
                float x_offset = x - (addedLeft ? 1.0 : 0.0);
                float y_offset = y - (addedFront ? 1.0 : 0.0);
                Vertex vertices[8] = {Vertex(x_offset  , y_offset  , z  , map[x  ][y  ][z  ]),
                                      Vertex(x_offset+1, y_offset  , z  , map[x+1][y  ][z  ]),
                                      Vertex(x_offset+1, y_offset+1, z  , map[x+1][y+1][z  ]),
                                      Vertex(x_offset  , y_offset+1, z  , map[x  ][y+1][z  ]),
                                      Vertex(x_offset  , y_offset  , z+1, map[x  ][y  ][z+1]),
                                      Vertex(x_offset+1, y_offset  , z+1, map[x+1][y  ][z+1]),
                                      Vertex(x_offset+1, y_offset+1, z+1, map[x+1][y+1][z+1]),
                                      Vertex(x_offset  , y_offset+1, z+1, map[x  ][y+1][z+1])
                                     };

                if (x == 0 && this->x == 0)
                {
                    vertices[0].isosurface = -1.0;
                    vertices[3].isosurface = -1.0;
                    vertices[4].isosurface = -1.0;
                    vertices[7].isosurface = -1.0;
                } else if (x == map.size() - (this->x > 0 ? 2 : 1) && this->lastChunkOnX)
                {
                    vertices[1].isosurface = -1.0;
                    vertices[2].isosurface = -1.0;
                    vertices[5].isosurface = -1.0;
                    vertices[6].isosurface = -1.0;
                }
                if (y == 0 && this->y == 0)
                {
                    vertices[0].isosurface = -1.0;
                    vertices[1].isosurface = -1.0;
                    vertices[4].isosurface = -1.0;
                    vertices[5].isosurface = -1.0;
                }if (y == map[x].size() - (this->y > 0 ? 2 : 1) && this->lastChunkOnY)
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
                    if (vertices[i].isosurface > isolevel)
                        cube_index ^= 1 << i;
                }
                int* edgesForTriangles = MarchingCubes::triangleTable[cube_index];
//                glPushName(reinterpret_cast<intptr_t>(grid->voxels[x][y][z]));
//                glBegin(GL_TRIANGLES);
                Vertex originalVertex;
                Vertex firstVertex;
                Vertex secondVertex;
                for (int i = 0; i < 16; i++) {
                    if (edgesForTriangles[i] == -1)
                        continue;
                    Vertex& v1 = vertices[MarchingCubes::edgeToCorner[edgesForTriangles[i]][0]];
                    Vertex& v2 = vertices[MarchingCubes::edgeToCorner[edgesForTriangles[i]][1]];

                    float interpolate = (isolevel - v1.isosurface) / (v2.isosurface - v1.isosurface);
                    Vertex midpoint = v1 - ((v1 - v2) * interpolate);
                    vertexArray.push_back(midpoint);
                    if (i % 3 == 0) {
                        originalVertex = midpoint;
                    }
                    else if (i % 3 == 1) {
                        firstVertex = midpoint;
                    }
                    else {
                        secondVertex = midpoint;
                        Vector3 normal = (firstVertex - originalVertex).cross((secondVertex - originalVertex)).normalize();
                    }
                }
            }
        }
    }
    return vertexArray;
}

void VoxelChunk::display(bool apply_marching_cubes, bool display_vertices, float isolevel)
{
    this->mesh.display();
    /*
    if (apply_marching_cubes) {
        MarchingCubes mc = MarchingCubes(*this);
        glPushMatrix();
        if (this->x == 0)
            glTranslatef(1.0, 0.0, 0.0);
        if (this->y == 0)
            glTranslatef(0.0, 1.0, 0.0);

        mc.display(isolevel);
        glPopMatrix();
    } else {
        for(int v_x = 0; v_x < sizeX; v_x++) {
            for(int v_y = 0; v_y < sizeY; v_y++) {
                for(int h = 0; h < height; h++) {
                    this->voxels[v_x][v_y][h]->display(apply_marching_cubes, display_vertices);
                }
            }
        }
    }*/
}

bool VoxelChunk::contains(Vector3 v) {
    return this->contains(v.x, v.y, v.z);
}

bool VoxelChunk::contains(float x, float y, float z) {
    return (this->x <= x && x < this->x + this->sizeX && this->y <= y && y < this->y + this->sizeY && 0 <= z && z < this->height);
}
