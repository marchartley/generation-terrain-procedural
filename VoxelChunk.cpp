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
VoxelChunk::~VoxelChunk()
{
    std::function<void(Voxel* v)> destruction([](Voxel* v) -> void {
        delete v;
    });
    this->applyToVoxels(destruction);
    this->voxels.clear();
}

void VoxelChunk::createMesh(bool applyMarchingCubes, bool updateMesh) {
    if (!needRemeshing)
        return;

    Voxel::currentLabelIndex = 1;
    Voxel::voxelGroups.clear();
    Voxel::voxelGroups.push_back(std::set<int>()); // First group reserved for the ones touching the ground
    this->computeGroups();

    if (!updateMesh)
        return;
    needRemeshing = false;

    std::vector<Vector3> colors;
    std::vector<Vector3> voxelVertices;
    if (applyMarchingCubes) {
        voxelVertices = this->applyMarchingCubes(true, &colors);
    }
    else {
        std::function<void(Voxel*)> createVoxels([&voxelVertices, &colors](Voxel* v) -> void {
            if ((bool)*v) {
                // Add the vertices to the global mesh
                std::vector<Vector3> vertice = v->getMeshVertices(true);
                voxelVertices.insert(voxelVertices.end(), vertice.begin(), vertice.end());
                // Add the colors to each vertex
                int X = 6; // Start with 6 faces
                for(auto& n : v->neighbors)
                    if (n.second && (bool)*n.second)
                        X--;    // Remove a face per neighbor
                X *= 6; // Multiply the number of face by the 6 vertex that defines it (2 triangles)
                for (int x = 0; x < X; x++) {
                    colors.push_back(Vector3((v->isOnGround ? 1.0 : 0.0), (v->isOnGround ? 0.0 : 1.0), 1.0));
//                        colors->push_back(HSVtoRGB((voxels[i][j][k]->group/((float)Voxel::voxelGroups.size()+1)), 1.0, 1.0));
                }
            }
        });
        this->applyToVoxels(createVoxels);
    }
    this->mesh.colorsArray = colors;
    this->mesh.fromArray(voxelVertices);
    this->mesh.update();
}

std::vector<Vector3> VoxelChunk::applyMarchingCubes(bool useGlobalCoords, std::vector<Vector3> *outColors)
{
    std::vector<Vector3> colors;
    std::vector<std::vector<std::vector<float> > > map = this->toFloat();
    /*for (int x = 0; x < this->sizeX; x++)
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
    }*/

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
                    colors.push_back(map[x][y][z] > 0.5 ? Vector3(150/255.0, 100/255.0, 30/255.0) : Vector3(224/255.0, 209/255.0, 72/255.0));
                    if (i % 3 == 0) {
                        originalVertex = midpoint;
                    }
                    else if (i % 3 == 1) {
                        firstVertex = midpoint;
                    }
                    else {
                        secondVertex = midpoint;
                        Vector3 mapOffset = (useGlobalCoords ? Vector3(this->x, this->y, 0.0) : Vector3(0.0, 0.0, 0.0));
                        vertexArray.push_back(originalVertex + mapOffset);
                        vertexArray.push_back(firstVertex + mapOffset);
                        vertexArray.push_back(secondVertex + mapOffset);
                    }
                }
            }
        }
    }
    if (outColors != nullptr)
        *outColors = colors;
    return vertexArray;
}

void VoxelChunk::makeItFall()
{
    this->applyToVoxels([](Voxel* v) -> void {
        Voxel* v_1 = v->neighbors[TOP];
        if (v_1 == nullptr) {
            v->isosurface = -0.01; // Just destroy the top voxels
            v->manual_isosurface = 0.0;
            return;
        }
        if (v->isOnGround || v_1->isOnGround)
            return;
        v->isosurface = v_1->isosurface;
        v->manual_isosurface = v_1->manual_isosurface;
    });
    this->needRemeshing = true;
    this->computeGroups();
}
void VoxelChunk::letGravityMakeSandFall()
{
    bool sandHadAction = false;
    this->applyToVoxels([&](Voxel* v) -> void {
        Voxel* v_1 = v->neighbors[TOP];
        if(v_1 && !*v_1) // If the top neighbor is air, don't care?
            return;
        if (v->getIsosurface() > 1.5)
            return;
        if(v->z == 0)
            return;
        if(v->getType() == DIRT) { // The bottom voxel is dirt
            return;
            if (v_1 == nullptr || v_1->getType() == DIRT)
                return;
            v->isosurface += v_1->isosurface;
            v->manual_isosurface += v_1->manual_isosurface;
            sandHadAction = true;
        } else { // The bottom voxel is sand or air
            if(!v_1 || v_1->getType() == DIRT) {
                // If top voxel is dirt or air,
                v->isosurface = 0.0;
                v->manual_isosurface = 0.0;
                sandHadAction = true;
                return;
            } else {
                v->isosurface += v_1->isosurface;
                v->manual_isosurface += v_1->manual_isosurface;
                sandHadAction = true;
            }
        }
        if (v_1) {
            v_1->isosurface = -0.01;
            v_1->manual_isosurface = -0.01;
        }
    });
    this->needRemeshing = sandHadAction;
//    this->computeGroups();
}

void VoxelChunk::resetVoxelsNeighbors() {
    this->applyToVoxels([](Voxel* v) -> void { v->resetNeighbors(); });
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
        if (v->getZ() == 0) {
            std::unordered_set<Voxel*> groundNeighbors;
            groundNeighbors.insert(v);
            while (groundNeighbors.size() != 0) {
                Voxel* n = (*groundNeighbors.begin());
                n->isOnGround = true;
                groundNeighbors.erase(groundNeighbors.begin());
                n->applyToNeighbors([&groundNeighbors](Voxel* v) -> void {
                    if(v != nullptr && !v->isOnGround && (bool)*v) {
                        groundNeighbors.insert(v);
                    }
                });
            }
        }
    });
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
std::vector<std::vector<std::vector<float>>>& VoxelChunk::toFloat() {
    this->voxelValues = std::vector<std::vector<std::vector<float>>>(this->sizeX);
    for(int x = 0; x < this->sizeX; x++) {
        this->voxelValues[x] = std::vector<std::vector<float>>(this->sizeY);
        for(int y = 0; y < this->sizeY; y++) {
            this->voxelValues[x][y] = std::vector<float>(this->height);
            for(int z = 0; z < this->height; z++) {
                this->voxelValues[x][y][z] = this->voxels[x][y][z]->getIsosurface();
            }
        }
    }
    return this->voxelValues;
}
std::vector<std::vector<std::vector<Voxel*>>>& VoxelChunk::toVoxels() {
    for(int x = 0; x < this->sizeX; x++) {
        for(int y = 0; y < this->sizeY; y++) {
            for(int z = 0; z < this->height; z++) {
                this->voxels[x][y][z]->isosurface = this->voxelValues[x][y][z];
            }
        }
    }
    return this->voxels;
}
