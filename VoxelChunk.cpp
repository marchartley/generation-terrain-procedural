#include "VoxelChunk.h"

#include <set>
#include <unordered_set>
#include <chrono>

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

VoxelChunk::VoxelChunk(int x, int y, int sizeX, int sizeY, int height, std::vector<std::vector<std::vector< float > > > iso_data, VoxelGrid* parent)
    : iso_data(iso_data), x(x), y(y), sizeX(sizeX), sizeY(sizeY), height(height), parent(parent) {

    this->voxels = std::vector<std::vector<std::vector<Voxel*>>>(this->sizeX);
    for(int v_x = 0; v_x < sizeX; v_x++) {
        this->voxels[v_x] = std::vector<std::vector<Voxel*>>(this->sizeY);
        for(int v_y = 0; v_y < sizeY; v_y++) {
            this->voxels[v_x][v_y] = std::vector<Voxel*>(this->height);
            for(int h = 0; h < height; h++) {
                Voxel* v = new Voxel(v_x, v_y, h, iso_data[v_x][v_y][h] > 0.0 ? DIRT : AIR, 1.0, iso_data[v_x][v_y][h], this);
                this->voxels[v_x][v_y][h] = v; //.push_back(v);
            }
        }
    }
}

VoxelChunk::VoxelChunk() : VoxelChunk(0, 0, 0, 0, 0, std::vector<std::vector<std::vector<float>>>(), nullptr)
{

}

void VoxelChunk::createMesh(bool applyMarchingCubes, bool updateMesh) {
    if (!needRemeshing)
        return;

    Voxel::currentLabelIndex = 1;
    Voxel::voxelGroups.clear();
    Voxel::voxelGroups.push_back(std::set<int>()); // First group reserved for the ones touching the ground
        /*this->applyToVoxels([](Voxel* v) -> void {
            v->group = -1;
            v->resetNeighbors();
        });*/
//    this->applyToVoxels([](Voxel* v) -> void {
//        if(!(bool)*v)
//            return;
//        v->group = v->getZ() == 0 ? 0 : -1; // If it's touching the ground, it's directly in first groupe
/*
        v->addNeighbor(v->getX() > 0 ? v->parent->voxels[v->getX()-1][v->getY()][v->getZ()] : nullptr);
        v->addNeighbor(v->getY() > 0 ? v->parent->voxels[v->getX()][v->getY()-1][v->getZ()] : nullptr);
        v->addNeighbor(v->getZ() > 0 ? v->parent->voxels[v->getX()][v->getY()][v->getZ()-1] : nullptr);

        if(v->getX() == 0 && v->parent->neighboring_chunks.find(LEFT) != v->parent->neighboring_chunks.end())
        {
            VoxelChunk* n = v->parent->neighboring_chunks[LEFT];
            v->addNeighbor(n->voxels[n->sizeX-1][v->getY()][v->getZ()]);
        }
        if(v->getY() == 0 && v->parent->neighboring_chunks.find(FRONT) != v->parent->neighboring_chunks.end())
        {
            VoxelChunk* n = v->parent->neighboring_chunks[FRONT];
            v->addNeighbor(n->voxels[v->getX()][n->sizeY - 1][v->getZ()]);
        }*/
/*
        if (v->group == -1)
        {
            v->group = Voxel::currentLabelIndex;
            Voxel::voxelGroups.push_back({Voxel::currentLabelIndex}); //std::vector<int>());
            Voxel::currentLabelIndex ++;
        } else {
            std::map<VOXEL_NEIGHBOR, Voxel*>::iterator it;
            for(it = v->neighbors.begin(); it != v->neighbors.end(); it++) {
                if(it->second && (bool)*it->second && it->second->group != -1)
                    Voxel::voxelGroups[v->group].insert(it->second->group);
            }
        }*/
//    });
    this->computeGroups();

    if (!updateMesh)
        return;
    needRemeshing = false;

    if (applyMarchingCubes) {
        this->mesh.fromArray(this->applyMarchingCubes());
    }
    else {
        std::vector<Vector3> voxelVertices;
        std::vector<Vector3> colors;
        this->applyToVoxels([voxelVertices=&voxelVertices, colors=&colors](Voxel* v) -> void {
            if ((bool)*v) {
                // Add the vertices to the global mesh
                std::vector<Vector3> vertice = v->getMeshVertices();
                voxelVertices->insert(voxelVertices->end(), vertice.begin(), vertice.end());
                // Add the colors to each vertex
                int X = 6; // Start with 6 faces
                for(std::map<VOXEL_NEIGHBOR, Voxel*>::iterator it = v->neighbors.begin(); it != v->neighbors.end(); it++)
                    if (it->second && (bool)*it->second)
                        X--;    // Remove a face per neighbor
                X *= 6; // Multiply the number of face by the 6 vertex that defines it (2 triangles)
                for (int x = 0; x < X; x++) {
                    colors->push_back(Vector3((v->isOnGround ? 1.0 : 0.0), (v->isOnGround ? 0.0 : 1.0), 1.0));
//                        colors->push_back(HSVtoRGB((voxels[i][j][k]->group/((float)Voxel::voxelGroups.size()+1)), 1.0, 1.0));
                }
            }
        });

        this->mesh.fromArray(voxelVertices);
        this->mesh.colorArrayFloat = Vector3::toArray(colors);
    }
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

    bool addedLeft = false;
    bool addedFront = false;
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

    if (addedLeft && addedFront) {
        VoxelChunk* n = this->neighboring_chunks[LEFT]->neighboring_chunks[FRONT];
        map[0].insert(map[0].begin(), std::vector<float>());
        for (int z = 0; z < this->height; z++)
            map[0][0].push_back(n->voxels[n->voxels.size() - 1][n->voxels[0].size() - 1][z]->getIsosurface());
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
                } else if (x == map.size() - 2 && this->lastChunkOnX) //(this->x > 0 ? 2 : 1) && this->lastChunkOnX)
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
                }if (y == map[x].size() - 2 && this->lastChunkOnY) // - (this->y > 0 ? 2 : 1) && this->lastChunkOnY)
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
//                        Vector3 normal = (firstVertex - originalVertex).cross((secondVertex - originalVertex)).normalize();
//                        normal.z *= -1.0;
//                        Voxel* in = this->parent->getVoxel(midpoint - normal);
//                        if(in)
//                            std::cout << in->getIsosurface() << " ";
                    }
                }
            }
        }
    }
    return vertexArray;
}

void VoxelChunk::makeItFall(int groupId)
{
    this->applyToVoxels([](Voxel* v) -> void {
        if (v->neighbors[TOP] == nullptr) {
            v->isosurface = -1.0; // Just destroy the top voxels
            v->manual_isosurface = 0.0;
            return;
        }
        Voxel* v_1 = v->neighbors[TOP];
        if (v->isOnGround || v_1->isOnGround)
            return;
        v->isosurface = v_1->isosurface;
        v->manual_isosurface = v_1->manual_isosurface;
    });
    this->needRemeshing = true;
    this->computeGroups();
}

void VoxelChunk::computeGroups()
{
    Voxel::currentLabelIndex = 1;
    Voxel::voxelGroups.clear();
    Voxel::voxelGroups.push_back(std::set<int>()); // First group reserved for the ones touching the ground

    this->applyToVoxels([](Voxel* v) -> void {
        v->group = -1;
        v->isOnGround = false;
    });

    this->applyToVoxels([](Voxel* v) -> void {
        if(!(bool)*v)
            return;
        v->group = v->getZ() == 0 ? 0 : -1; // If it's touching the ground, it's directly in first groupe
        /*
        v->addNeighbor(v->getX() > 0 ? v->parent->voxels[v->getX()-1][v->getY()][v->getZ()] : nullptr);
        v->addNeighbor(v->getY() > 0 ? v->parent->voxels[v->getX()][v->getY()-1][v->getZ()] : nullptr);
        v->addNeighbor(v->getZ() > 0 ? v->parent->voxels[v->getX()][v->getY()][v->getZ()-1] : nullptr);

        if(v->getX() == 0 && v->parent->neighboring_chunks.find(LEFT) != v->parent->neighboring_chunks.end())
        {
            VoxelChunk* n = v->parent->neighboring_chunks[LEFT];
            v->addNeighbor(n->voxels[n->sizeX-1][v->getY()][v->getZ()]);
        }
        if(v->getY() == 0 && v->parent->neighboring_chunks.find(FRONT) != v->parent->neighboring_chunks.end())
        {
            VoxelChunk* n = v->parent->neighboring_chunks[FRONT];
            v->addNeighbor(n->voxels[v->getX()][n->sizeY - 1][v->getZ()]);
        }*//*
        if (v->group == -1)
        {
            v->group = Voxel::currentLabelIndex;
            Voxel::voxelGroups.push_back({Voxel::currentLabelIndex}); //std::vector<int>());
            Voxel::currentLabelIndex ++;
        } else {
            std::map<VOXEL_NEIGHBOR, Voxel*>::iterator it;
            for(it = v->neighbors.begin(); it != v->neighbors.end(); it++) {
                if(it->second && it->second->getType()!=TerrainTypes::AIR && it->second->group != -1)
                    Voxel::voxelGroups[v->group].insert(it->second->group);
            }
        }*/
        if (v->getZ() == 0) {
            std::unordered_set<Voxel*> groundNeighbors;
            groundNeighbors.insert(v);
            while (groundNeighbors.size() != 0) {
                Voxel* n = (*groundNeighbors.begin());
                n->isOnGround = true;
                groundNeighbors.erase(groundNeighbors.begin());
                for(std::map<VOXEL_NEIGHBOR, Voxel*>::iterator it = n->neighbors.begin(); it != n->neighbors.end(); it++)
                    if(it->second && (bool)*it->second && !it->second->isOnGround)
                        groundNeighbors.insert(it->second);
            }
        }
    });
/*
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
    }*/
}


void VoxelChunk::display()
{
    this->mesh.display();
}

bool VoxelChunk::contains(Vector3 v) {
    return this->contains(v.x, v.y, v.z);
}

bool VoxelChunk::contains(float x, float y, float z) {
    return (this->x <= x && x < this->x + this->sizeX && this->y <= y && y < this->y + this->sizeY && 0 <= z && z < this->height);
}
