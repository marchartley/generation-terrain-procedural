#include "VoxelChunk.h"

#include <set>
Vector3 HSVtoRGB(float H, float S,float V){
    H *= 360;
    float s = S;
    float v = V;
    float C = s*v;
    float X = C*(1-abs(fmod(H/60.0, 2)-1));
    float m = v-C;
    float r,g,b;
    if(H >= 0 && H < 60){
        r = C,g = X,b = 0;
    }
    else if(H >= 60 && H < 120){
        r = X,g = C,b = 0;
    }
    else if(H >= 120 && H < 180){
        r = 0,g = C,b = X;
    }
    else if(H >= 180 && H < 240){
        r = 0,g = X,b = C;
    }
    else if(H >= 240 && H < 300){
        r = X,g = 0,b = C;
    }
    else{
        r = C,g = 0,b = X;
    }
    float R = (r+m);
    float G = (g+m);
    float B = (b+m);
//    std::cout << Vector3(H, S, V) << " " << Vector3(R, G, B) << std::endl;
    return Vector3(R, G, B);
}

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

void VoxelChunk::createMesh(bool updateMesh) {
    if (!needRemeshing)
        return;
    Voxel::currentLabelIndex = 1;
    Voxel::voxelGroups.clear();
    Voxel::voxelGroups.push_back(std::set<int>()); // First group reserved for the ones touching the ground
    if (this->voxels.size() == 0) {
        this->voxels = std::vector<std::vector<std::vector<Voxel*>>>();
        for(int v_x = 0; v_x < sizeX; v_x++) {
            this->voxels.push_back(std::vector<std::vector<Voxel*>>());
            for(int v_y = 0; v_y < sizeY; v_y++) {
                this->voxels[v_x].push_back(std::vector<Voxel*>());
                for(int h = 0; h < height; h++) {
                    Voxel* v = new Voxel(v_x, v_y, h, this->data[v_x][v_y][h], 1.0, iso_data[v_x][v_y][h]);
                    v->parent = this;
                    this->voxels[v_x][v_y].push_back(v);
                }
            }
        }
    } else {
        for(int v_x = 0; v_x < sizeX; v_x++) {
            for(int v_y = 0; v_y < sizeY; v_y++) {
                for(int h = 0; h < height; h++) {
                    Voxel* v = this->voxels[v_x][v_y][h];
                    v->group = -1;
                    v->resetNeighbors();
                }
            }
        }
    }

    for(int v_x = 0; v_x < sizeX; v_x++) {
        for(int v_y = 0; v_y < sizeY; v_y++) {
            for(int h = 0; h < height; h++) {
                Voxel* v = this->voxels[v_x][v_y][h];
//                v->resetNeighbors();
                v->type = (v->getIsosurface() > 0.0 ? TerrainTypes::DIRT : TerrainTypes::AIR);
                if(v->type == TerrainTypes::AIR)
                    continue;
                v->group = h == 0 ? 0 : -1; // If it's touching the ground, it's directly in first groupe
                v->addNeighbor(v_x > 0 ? this->voxels[v_x-1][v_y][h] : nullptr);
                v->addNeighbor(v_y > 0 ? this->voxels[v_x][v_y-1][h] : nullptr);
                v->addNeighbor(h > 0 ? this->voxels[v_x][v_y][h-1] : nullptr);

                if(v_x == 0 && this->neighboring_chunks.find(LEFT) != this->neighboring_chunks.end())
                {
                    VoxelChunk* n = this->neighboring_chunks[LEFT];
                    v->addNeighbor(n->voxels[n->sizeX-1][v_y][h]);
                }
                if(v_y == 0 && this->neighboring_chunks.find(FRONT) != this->neighboring_chunks.end())
                {
                    VoxelChunk* n = this->neighboring_chunks[FRONT];
                    v->addNeighbor(n->voxels[v_x][n->sizeY - 1][h]);
                }

                if (v->group == -1)
                {
                    v->group = Voxel::currentLabelIndex;
                    Voxel::voxelGroups.push_back({Voxel::currentLabelIndex}); //std::vector<int>());
                    Voxel::currentLabelIndex ++;
                } else {
                    std::map<VOXEL_NEIGHBOR, Voxel*>::iterator it;
                    for(it = v->neighbors.begin(); it != v->neighbors.end(); it++) {
                        if(it->second && it->second->type!=TerrainTypes::AIR && it->second->group != -1)
                            Voxel::voxelGroups[v->group].insert(it->second->group);
                    }
                }
            }
        }
    }
    this->computeGroups();

    if (!updateMesh)
        return;
    needRemeshing = false;

    std::vector<Vector3> voxelVertices;
    for (int i = 0; i < voxels.size(); i++)
        for (int j = 0; j < voxels[i].size(); j++)
            for (int k = 0; k < voxels[i][j].size(); k++) {
                std::vector<Vector3> vertice = voxels[i][j][k]->getMeshVertices();
                voxelVertices.insert(voxelVertices.end(), vertice.begin(), vertice.end());
            }
    this->mesh.fromArray(voxelVertices);

    std::vector<Vector3> colors;
    for (int i = 0; i < voxels.size(); i++)
        for (int j = 0; j < voxels[i].size(); j++)
            for (int k = 0; k < voxels[i][j].size(); k++) {
                Voxel* v = voxels[i][j][k];
                if (voxels[i][j][k]->type != TerrainTypes::AIR) {
                    int X = 6;
                    for(std::map<VOXEL_NEIGHBOR, Voxel*>::iterator it = v->neighbors.begin(); it != v->neighbors.end(); it++)
                        if (it->second)
                            X--;
//                    X -= (i == 0 ? 1:0) + (j == 0 ? 1:0) + (k == 0 ? 1:0) + (i == voxels.size()-1 ? 1:0) + (j == voxels.size()-1 ? 1:0) + (k == voxels.size()-1 ? 1:0);
                    X *= 6;
                    for (int x = 0; x < X; x++) {
                        colors.push_back(HSVtoRGB((voxels[i][j][k]->group/((float)Voxel::voxelGroups.size()+1)), 1.0, 1.0));
                    }
                }
            }
    this->mesh.randomStuff = Vector3::toArray(colors);
//    std::cout << this->mesh.normalsArrayFloat.size() << " -- " << this->mesh.randomStuff.size() << std::endl;

//    this->mesh.fromArray(this->applyMarchingCubes());
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

void VoxelChunk::makeItFall(int groupId)
{
//    if (groupId == -1) {
//        // Don't do it for groupe 1
//        for (size_t i = 1; i < Voxel::voxelGroups.size(); i++)
//            this->makeItFall(i);
//        return;
//    }
    for (int i = 0; i < voxels.size(); i++) {
        for (int j = 0; j < voxels[i].size(); j++) {
            for (int k = 0; k < voxels[i][j].size()-1; k++) {
                Voxel* v_1 = voxels[i][j][k+1];
                Voxel* v = voxels[i][j][k];
//                if (v->type == TerrainTypes::AIR)
//                    continue;
//                if(v->group != groupId)
//                    continue;
//                if(v->has_neighbors[BOTTOM] || k == 0)
//                    continue;
//                if (k == 0)
//                    continue;
                if (v->group == 0 || v_1->group == 0)
                    continue;
                if (k+1 == this->height-1)
                    v_1->isosurface = -1.0;
                v->isosurface = v_1->isosurface;
                /*
                float tmpIso = voxels[i][j][k]->isosurface, tmpManualIso = voxels[i][j][k+1]->manual_isosurface;
                int tmpGroup = voxels[i][j][k]->group;
                voxels[i][j][k]->isosurface = v->isosurface;
                voxels[i][j][k]->manual_isosurface = v->manual_isosurface;
                voxels[i][j][k]->group = v->group;
                if (k == this->height) {
                    v->isosurface = -1000;
                    v->manual_isosurface = -1000;
                    v->group = -1;
                } else {
                    v->isosurface = tmpIso;// * .1;
                    v->manual_isosurface = tmpManualIso;
                    v->group = tmpGroup;
                }
                v->resetNeighbors();*/
            }
        }
    }
//    for (int i = 0; i < voxels.size(); i++) {
//        for (int j = 0; j < voxels[i].size(); j++) {
//            for (int k = 0; k < voxels[i][j].size(); k++) {
//                voxels[i][j][k]->resetNeighbors();
//            }
//        }
//    }
    this->needRemeshing = true;
    this->computeGroups();
}

void VoxelChunk::computeGroups()
{
    Voxel::currentLabelIndex = 1;
    Voxel::voxelGroups.clear();
    Voxel::voxelGroups.push_back(std::set<int>()); // First group reserved for the ones touching the ground
    for(int v_x = 0; v_x < sizeX; v_x++) {
        for(int v_y = 0; v_y < sizeY; v_y++) {
            for(int h = 0; h < height; h++) {
                Voxel* v = this->voxels[v_x][v_y][h];
                v->group = -1;
                v->resetNeighbors();
                if (v->type == AIR)
                    continue;
                v->group = h == 0 ? 0 : -1; // If it's touching the ground, it's directly in first groupe
                v->addNeighbor(v_x > 0 ? this->voxels[v_x-1][v_y][h] : nullptr);
                v->addNeighbor(v_y > 0 ? this->voxels[v_x][v_y-1][h] : nullptr);
                v->addNeighbor(h > 0 ? this->voxels[v_x][v_y][h-1] : nullptr);

                if(v_x == 0 && this->neighboring_chunks.find(LEFT) != this->neighboring_chunks.end())
                {
                    VoxelChunk* n = this->neighboring_chunks[LEFT];
                    v->addNeighbor(n->voxels[n->sizeX-1][v_y][h]);
                }
                if(v_y == 0 && this->neighboring_chunks.find(FRONT) != this->neighboring_chunks.end())
                {
                    VoxelChunk* n = this->neighboring_chunks[FRONT];
                    v->addNeighbor(n->voxels[v_x][n->sizeY - 1][h]);
                }
                if (v->group == -1)
                {
                    v->group = Voxel::currentLabelIndex;
                    Voxel::voxelGroups.push_back({Voxel::currentLabelIndex}); //std::vector<int>());
                    Voxel::currentLabelIndex ++;
                } else {
                    std::map<VOXEL_NEIGHBOR, Voxel*>::iterator it;
                    for(it = v->neighbors.begin(); it != v->neighbors.end(); it++) {
                        if(it->second && it->second->type!=TerrainTypes::AIR && it->second->group != -1)
                            Voxel::voxelGroups[v->group].insert(it->second->group);
                    }
                }
            }
        }
    }

    std::vector<std::set<int>> tempVoxelGroup;
    for(int i = 0; i < Voxel::voxelGroups.size(); i++)
    {
        int found = -1;
        for (std::set<int>::iterator i_lab = Voxel::voxelGroups[i].begin(); i_lab != Voxel::voxelGroups[i].end(); i_lab++)
        {
            for(int ii = 0; ii < tempVoxelGroup.size(); ii++)
            {
                if(tempVoxelGroup[ii].find(*i_lab) != tempVoxelGroup[ii].end()) {
                    found = ii;
                    break;
                }
            }
        }
        if (found == -1) {
            tempVoxelGroup.push_back(Voxel::voxelGroups[i]);
        } else {
            tempVoxelGroup[found].insert(Voxel::voxelGroups[i].begin(), Voxel::voxelGroups[i].end());
        }
    }
    Voxel::voxelGroups = tempVoxelGroup;

    for (int i = 0; i < voxels.size(); i++) {
        for (int j = 0; j < voxels[i].size(); j++) {
            for (int k = 0; k < voxels[i][j].size(); k++) {
                Voxel* v = voxels[i][j][k];
                int _ii = 0;
                for(_ii = 0; _ii < Voxel::voxelGroups.size(); _ii++) {
                    if(Voxel::voxelGroups[_ii].find(v->group) != Voxel::voxelGroups[_ii].end()) {
                        v->group = _ii;
                        break;
                    }
                }
            }
        }
    }
}


void VoxelChunk::display(bool apply_marching_cubes, bool display_vertices, float isolevel)
{
    this->mesh.display();
}

bool VoxelChunk::contains(Vector3 v) {
    return this->contains(v.x, v.y, v.z);
}

bool VoxelChunk::contains(float x, float y, float z) {
    return (this->x <= x && x < this->x + this->sizeX && this->y <= y && y < this->y + this->sizeY && 0 <= z && z < this->height);
}
