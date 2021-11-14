#include "VoxelChunk.h"

VoxelChunk::VoxelChunk(int _x, int _y, int _sizeX, int _sizeY, int _height, std::vector<std::vector<std::vector< TerrainTypes > > > _data)
    : data(_data), height(_height), x(_x), y(_y), sizeX(_sizeX), sizeY(_sizeY) {
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
VoxelChunk::VoxelChunk(int x, int y, int sizeX, int sizeY, int height, std::vector<std::vector<std::vector< float > > > iso_data)
    : iso_data(iso_data), height(height), x(x), y(y), sizeX(sizeX), sizeY(sizeY) {
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

VoxelChunk::VoxelChunk() : VoxelChunk(0, 0, 0, 0, 0, std::vector<std::vector<std::vector<TerrainTypes>>>())
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

//                    if(v_x > 0) {
//                        v->isosurfaces[0] = this->voxels[v_x-1][v_y][h]->isosurfaces[1];
//                        v->isosurfaces[3] = this->voxels[v_x-1][v_y][h]->isosurfaces[2];
//                        v->isosurfaces[4] = this->voxels[v_x-1][v_y][h]->isosurfaces[5];
//                        v->isosurfaces[7] = this->voxels[v_x-1][v_y][h]->isosurfaces[6];
//                    }
//                    if(v_y > 0) {
//                        v->isosurfaces[2] = this->voxels[v_x][v_y-1][h]->isosurfaces[0];
//                        v->isosurfaces[3] = this->voxels[v_x][v_y-1][h]->isosurfaces[1];
//                        v->isosurfaces[6] = this->voxels[v_x][v_y-1][h]->isosurfaces[5];
//                        v->isosurfaces[7] = this->voxels[v_x][v_y-1][h]->isosurfaces[4];
//                    }
//                    if(h > 0) {
//                        v->isosurfaces[0] = this->voxels[v_x][v_y][h-1]->isosurfaces[4];
//                        v->isosurfaces[1] = this->voxels[v_x][v_y][h-1]->isosurfaces[5];
//                        v->isosurfaces[2] = this->voxels[v_x][v_y][h-1]->isosurfaces[6];
//                        v->isosurfaces[3] = this->voxels[v_x][v_y][h-1]->isosurfaces[7];
//                    }

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
//                }
            }
        }
    }
    needRemeshing = false;





    /*for(int v_x = 0; v_x < sizeX; v_x++) {
        for(int v_y = 0; v_y < sizeY; v_y++) {
            for(int h = 0; h < height; h++) {
                Voxel* v = this->voxels[v_x][v_y][h];
                if (*v || !*v) {
                    if(v == NULL)
                        continue;
//                    for(int i = 0; i < 8; i++) {
//                        v->vertices[i].isosurface = *v ? 1.0 : -1.0;
//                    }
                    // Compute isosurfaces
                    /*if (v_x == 0){
                        v->vertices[0].isosurface = -1.0;
                        v->vertices[3].isosurface = -1.0;
                        v->vertices[4].isosurface = -1.0;
                        v->vertices[7].isosurface = -1.0;
                    } else if (v_x == sizeX - 1) {
                        v->vertices[1].isosurface = -1.0;
                        v->vertices[2].isosurface = -1.0;
                        v->vertices[5].isosurface = -1.0;
                        v->vertices[6].isosurface = -1.0;
                    }

                    if (v_y == 0){
                        v->vertices[2].isosurface = -1.0;
                        v->vertices[3].isosurface = -1.0;
                        v->vertices[6].isosurface = -1.0;
                        v->vertices[7].isosurface = -1.0;
                    } else if (v_y == sizeY - 1) {
                        v->vertices[0].isosurface = -1.0;
                        v->vertices[1].isosurface = -1.0;
                        v->vertices[4].isosurface = -1.0;
                        v->vertices[5].isosurface = -1.0;
                    }

                    if (h == 0){
                        v->vertices[0].isosurface = -1.0;
                        v->vertices[1].isosurface = -1.0;
                        v->vertices[2].isosurface = -1.0;
                        v->vertices[3].isosurface = -1.0;
                    } else if (h == height - 1) {
                        v->vertices[4].isosurface = -1.0;
                        v->vertices[5].isosurface = -1.0;
                        v->vertices[6].isosurface = -1.0;
                        v->vertices[7].isosurface = -1.0;
                    }*//*
                    if (v_x > 0) {
                        if (v_y > 0) {
                            if (h > 0) {
                                v->vertices[3].isosurface = (*this->voxels[v_x - 1][v_y - 1][h - 1] ? 1.0 : -1.0);
                            }
                            if (h < height - 1) {
                                v->vertices[7].isosurface = (*this->voxels[v_x - 1][v_y - 1][h + 1] ? 1.0 : -1.0);
                            }
                        }
                        if (v_y < sizeY - 1) {
                            if (h > 0) {
                                v->vertices[0].isosurface = (*this->voxels[v_x - 1][v_y + 1][h - 1] ? 1.0 : -1.0);
                            }
                            if (h < height - 1) {
                                v->vertices[4].isosurface = (*this->voxels[v_x - 1][v_y + 1][h + 1] ? 1.0 : -1.0);
                            }
                        }
                    }
                    if (v_x < sizeX - 1) {
                        if (v_y > 0) {
                            if (h > 0) {
                                v->vertices[2].isosurface = (*this->voxels[v_x + 1][v_y - 1][h - 1] ? 1.0 : -1.0);
                            }
                            if (h < height - 1) {
                                v->vertices[6].isosurface = (*this->voxels[v_x + 1][v_y - 1][h + 1] ? 1.0 : -1.0);
                            }
                        }
                        if (v_y < sizeY - 1) {
                            if (h > 0) {
                                v->vertices[1].isosurface = (*this->voxels[v_x + 1][v_y + 1][h - 1] ? 1.0 : -1.0);
                            }
                            if (h < height - 1) {
                                v->vertices[5].isosurface = (*this->voxels[v_x + 1][v_y + 1][h + 1] ? 1.0 : -1.0);
                            }
                        }
                    }*//*
                    int a = 0,
                        b = 1,
                        c = 2,
                        d = 3;
                    *v->isosurfaces[a] = (v->has_neighbors[LEFT] && v->has_neighbors[BOTTOM] && v->has_neighbors[BACK]) ? 1.0 : -1.0;
                    *v->isosurfaces[b] = (v->has_neighbors[RIGHT] && v->has_neighbors[BOTTOM] && v->has_neighbors[BACK]) ? 1.0 : -1.0;
                    *v->isosurfaces[c] = (v->has_neighbors[RIGHT] && v->has_neighbors[BOTTOM] && v->has_neighbors[FRONT]) ? 1.0 : -1.0;
                    *v->isosurfaces[d] = (v->has_neighbors[LEFT] && v->has_neighbors[BOTTOM] && v->has_neighbors[FRONT]) ? 1.0 : -1.0;
                    *v->isosurfaces[a+4] = (v->has_neighbors[LEFT] && v->has_neighbors[TOP] && v->has_neighbors[BACK]) ? 1.0 : -1.0;
                    *v->isosurfaces[b+4] = (v->has_neighbors[RIGHT] && v->has_neighbors[TOP] && v->has_neighbors[BACK]) ? 1.0 : -1.0;
                    *v->isosurfaces[c+4] = (v->has_neighbors[RIGHT] && v->has_neighbors[TOP] && v->has_neighbors[FRONT]) ? 1.0 : -1.0;
                    *v->isosurfaces[d+4] = (v->has_neighbors[LEFT] && v->has_neighbors[TOP] && v->has_neighbors[FRONT]) ? 1.0 : -1.0;
                }
            }
        }
    }*/
}
void VoxelChunk::display(bool apply_marching_cubes, bool display_vertices, float isolevel)
{
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
    }
}

bool VoxelChunk::contains(Vector3 v) {
    return this->contains(v.x, v.y, v.z);
}

bool VoxelChunk::contains(float x, float y, float z) {
    return (this->x <= x && x < this->x + this->sizeX && this->y <= y && y < this->y + this->sizeY && 0 <= z && z < this->height);
}
