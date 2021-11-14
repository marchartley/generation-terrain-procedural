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

    MarchingCubes mc = MarchingCubes(*this);
    mc.createMesh();
    this->mesh = mc.mesh;

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
