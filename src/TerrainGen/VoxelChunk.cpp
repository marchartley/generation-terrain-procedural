#include "TerrainGen/VoxelChunk.h"

#include <set>
#include <unordered_set>
#include <chrono>
#include "Utils/Utils.h"

VoxelChunk::VoxelChunk(int x, int y, int sizeX, int sizeY, int height, Matrix3<float> iso_data, std::shared_ptr<VoxelGrid> parent)
    : iso_data(iso_data), x(x), y(y), sizeX(sizeX), sizeY(sizeY), sizeZ(height), parent(parent) {

    this->voxelGroups = Matrix3<int>(this->sizeX, this->sizeY, this->sizeZ, -1);
    this->applyModification(iso_data);
    this->flowField = Matrix3<Vector3>(this->sizeX, this->sizeY, this->sizeZ);
    this->updateLoDsAvailable();
    this->computeDistanceField();
}

VoxelChunk::VoxelChunk() : VoxelChunk(0, 0, 0, 0, 0, Matrix3<float>(), nullptr)
{

}
VoxelChunk::~VoxelChunk()
{
}

void VoxelChunk::computeFlowfield(Vector3 sea_current)
{
    computeDistanceField();
    Matrix3<Vector3> pressionGradient = this->pressureField.gradient();
    this->flowField = pressionGradient * (-1.f) + sea_current;
}

void VoxelChunk::computeDistanceField()
{
    this->distanceField = Matrix3<int>(this->sizeX, this->sizeY, this->sizeZ, 9999999);
    this->pressureField = Matrix3<float>(this->sizeX, this->sizeY, this->sizeZ, 0);
    // Check left
    for (int x = 1; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                if (this->getVoxelValue(x, y, z) > 0) {
                    distanceField(x, y, z) = 0;
                } else {
                    distanceField(x, y, z) = std::min(distanceField(x, y, z), distanceField(x - 1, y, z) + 1);
                }
            }
        }
    }
    // Check right
    for (int x = 0; x < this->sizeX - 1; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                if (this->getVoxelValue(x, y, z) > 0) {
                    distanceField(x, y, z) = 0;
                } else {
                    distanceField(x, y, z) = std::min(distanceField(x, y, z), distanceField(x + 1, y, z) + 1);
                }
            }
        }
    }
    // Check front
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 1; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                if (this->getVoxelValue(x, y, z) > 0) {
                    distanceField(x, y, z) = 0;
                } else {
                    distanceField(x, y, z) = std::min(distanceField(x, y, z), distanceField(x, y - 1, z) + 1);
                }
            }
        }
    }
    // Check back
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY - 1; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                if (this->getVoxelValue(x, y, z) > 0) {
                    distanceField(x, y, z) = 0;
                } else {
                    distanceField(x, y, z) = std::min(distanceField(x, y, z), distanceField(x, y + 1, z) + 1);
                }
            }
        }
    }
    // Check down
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 1; z < this->sizeZ; z++) {
                if (this->getVoxelValue(x, y, z) > 0) {
                    distanceField(x, y, z) = 0;
                } else {
                    distanceField(x, y, z) = std::min(distanceField(x, y, z), distanceField(x, y, z - 1) + 1);
                }
            }
        }
    }
    // Check up
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ - 1; z++) {
                if (this->getVoxelValue(x, y, z) > 0) {
                    distanceField(x, y, z) = 0;
                } else {
                    distanceField(x, y, z) = std::min(distanceField(x, y, z), distanceField(x, y, z + 1) + 1);
                }
            }
        }
    }
    // compute pressure ( 1 / distance^2 ?)
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                if (distanceField(x, y, z) > 0) // Works better with 1/dist than 1/dist^2
                    this->pressureField(x, y, z) = 1 / (float)this->distanceField(x, y, z); //std::min(1.f , 1 / ((float)std::pow(this->distanceField(x, y, z), 2)));
                else
                    this->pressureField(x, y, z) = 1.0;
            }
        }
    }
}

void VoxelChunk::updateLoDsAvailable() {
    this->LoDs.clear();
    int x_add = 1; //(this->x == 0 ? 0 : 0);
    int y_add = 1; //(this->y == 0 ? 0 : 0);
    int z_add = 1;
    for (int i = 1; i < std::min(std::min(this->sizeX, this->sizeY), this->sizeZ); i++)
        if ((this->sizeX-x_add) % i == 0 && (this->sizeY-y_add) % i == 0 && (this->sizeZ-z_add) % i == 0)
            this->LoDs.push_back(i);
}
void VoxelChunk::createMesh(bool applyMarchingCubes, bool updateMesh) {
    if (!needRemeshing)
        return;

    if (!updateMesh)
        return;
    needRemeshing = false;

    std::vector<Vector3> colors;
    std::vector<Vector3> normals;
    std::vector<Vector3> voxelVertices;
    if (applyMarchingCubes) {
        voxelVertices = this->applyMarchingCubes(true, &colors, &normals);
        this->mesh.colorsArray = colors;
        this->mesh.normalsArray = normals;
        this->mesh.fromArray(voxelVertices);
        this->mesh.update();
    }
    else {/*
        this->mesh.fromArray(CubeMesh::cubesVertices);
        std::vector<Vector3> voxelsPos;
        Matrix3<float> allVals = this->getVoxelValues();
        for (int x = 0; x < this->sizeX; x++) {
            for (int y = 0; y < this->sizeY; y++) {
                for (int z = 0; z < this->height; z++) {
                    if (allVals.at(x, y, z) > 0.f)
                        voxelsPos.push_back(Vector3(x, y, z));
                }
            }
        }
        this->mesh.update(voxelsPos);*/
    }
}

Vector3 VoxelChunk::computeNormal(Vector3 pos)
{
    return this->computeNormal(pos.x, pos.y, pos.z);
}
Vector3 VoxelChunk::computeNormal(int x, int y, int z)
{
    Vertex vertices[8] = {Vertex(x    , y    , z    , this->parent->getVoxelValue(x    , y    , z    )),
                          Vertex(x + 1, y    , z    , this->parent->getVoxelValue(x + 1, y    , z    )),
                          Vertex(x + 1, y + 1, z    , this->parent->getVoxelValue(x + 1, y + 1, z    )),
                          Vertex(x    , y + 1, z    , this->parent->getVoxelValue(x    , y + 1, z    )),
                          Vertex(x    , y    , z + 1, this->parent->getVoxelValue(x    , y    , z + 1)),
                          Vertex(x + 1, y    , z + 1, this->parent->getVoxelValue(x + 1, y    , z + 1)),
                          Vertex(x + 1, y + 1, z + 1, this->parent->getVoxelValue(x + 1, y + 1, z + 1)),
                          Vertex(x    , y + 1, z + 1, this->parent->getVoxelValue(x    , y + 1, z + 1))
                         };

    std::vector<Vector3> tris = this->computeMarchingCube(vertices, 0.0, false);
    Vector3 normal;
    for(size_t i = 0; i < tris.size(); i += 3)
    {
        normal += (tris[i+2] - tris[i]).cross((tris[i+1] - tris[i])) +
                    (tris[i] - tris[i+1]).cross((tris[i+2] - tris[i+1])) +
                    (tris[i+1] - tris[i+2]).cross((tris[i] - tris[i+2]));
    }
    if (tris.size() == 0)
        return normal;
    if (this->parent->getVoxelValue(Vector3(x, y, z) + normal) > 0)
        normal *= -1.f;
    return normal.normalize();

}

std::vector<Vector3> VoxelChunk::computeMarchingCube(Vertex vertices[8], float isolevel, bool useGlobalCoords, std::vector<Vector3>* outColors, std::vector<Vector3>* outNormals)
{
    std::vector<Vector3> vertexArray;
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
        if (outNormals != nullptr) {
            outNormals->push_back(Vector3((v2 - v1) * (v1.isosurface > v2.isosurface ? 1.f : -1.f))); //this->computeNormal(vertices[0]));
        }
        if (outColors != nullptr) {
//            float blueVal = this->parent->getVoxelValue(v1) - this->parent->getOriginalVoxelValue(v1) + this->parent->getVoxelValue(v2) - this->parent->getOriginalVoxelValue(v2);
//            blueVal = abs(blueVal / 2.0);$
            float blueVal = 30/255.0;
            outColors->push_back(vertices[0].isosurface > 0.5 ?  //map(x, y, z) > 0.5 ?
                                 Vector3(150/255.0, 100/255.0, blueVal) : //30/255.0) :
                                 Vector3(224/255.0, 209/255.0, blueVal)); //72/255.0));
        }
        if (i % 3 == 0) {
            originalVertex = midpoint;
        }
        else if (i % 3 == 1) {
            firstVertex = midpoint;
        }
        else {
            secondVertex = midpoint;
            Vector3 mapOffset = (useGlobalCoords ? Vector3(this->x, this->y, 0.0) : Vector3(0.0, 0.0, 0.0));
            vertexArray.insert(vertexArray.end(), {originalVertex + mapOffset,
                                                   firstVertex + mapOffset,
                                                   secondVertex + mapOffset});
        }
    }
    return vertexArray;
}
std::vector<Vector3> VoxelChunk::applyMarchingCubes(bool useGlobalCoords, std::vector<Vector3>* outColors, std::vector<Vector3>* outNormals)
{
    int LoD = this->LoDs[this->LoDIndex % this->LoDs.size()];//std::max(0, std::min(LoD, std::min(std::min(this->sizeX, this->sizeY), this->height)));
    std::vector<Vector3> colors;
    std::vector<Vector3> normals;
    Matrix3<float> map = this->getVoxelValues();

//    map = map.meanSmooth(3, 3, 3, true);

    bool addedLeft = false;
    bool addedFront = false;
    if (this->neighboring_chunks.find(LEFT) != this->neighboring_chunks.end()) {
        std::shared_ptr<VoxelChunk> n = this->neighboring_chunks[LEFT];
        map.insertRow(0, 0); // Insert column at the begining of the X-axis
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++)
            {
                map(0, y, z) = n->getVoxelValue(n->sizeX - 1, y, z);
            }
        }
        addedLeft = true;
    }

    if (this->neighboring_chunks.find(FRONT) != this->neighboring_chunks.end()) {
        std::shared_ptr<VoxelChunk> n = this->neighboring_chunks[FRONT];
        int offset = addedLeft ? 1 : 0;
        map.insertRow(/*offset*/0, 1);
        for (int x = 0; x < this->sizeX; x++) {
            for (int z = 0; z < this->sizeZ; z++)
            {
                map(x + offset, 0, z) = n->getVoxelValue(x, n->sizeY - 1, z);
            }
        }
        addedFront = true;
    }
    if (addedLeft && addedFront) {
        std::shared_ptr<VoxelChunk> n = this->neighboring_chunks[LEFT]->neighboring_chunks[FRONT];
        for (int z = 0; z < this->sizeZ; z++)
            map(0, 0, z) = n->getVoxelValue(n->sizeX - 1, n->sizeY - 1, z);
    }


    if (this->x == 0) {
        for (int y = 0; y < map.sizeY; y++)
            for (int z = 0; z < this->sizeZ; z++)
                map(0, y, z) = -1;
    }
    if (this->lastChunkOnX)
        for (int y = 0; y < map.sizeY; y++)
            for (int z = 0; z < this->sizeZ; z++)
                map(map.sizeX - 1, y, z) = -1;
    if (this->y == 0) {
        for (int x = 0; x < map.sizeX; x++)
            for (int z = 0; z < this->sizeZ; z++)
                map(x, 0, z) = -1;
    }
    if (this->lastChunkOnY)
        for (int x = 0; x < map.sizeX; x++)
            for (int z = 0; z < this->sizeZ; z++)
                map(x, map.sizeY - 1, z) = -1;

    for (int x = 0; x < map.sizeX; x++) {
        for (int y = 0; y < map.sizeY; y++) {
            map(x, y, 0) = -1;
            map(x, y, map.sizeZ - 1) = -1;
        }
    }
    float isolevel = 0.0;
    std::vector<Vector3> vertexArray;
    int x = 0, y = 0, z = 0;
    for (x = 0; x < map.sizeX - 1; x += LoD) {
        for (y = 0; y < map.sizeY - 1; y += LoD) {
            for (z = 0; z < map.sizeZ - 1; z += LoD) {
                int LoDX = std::min(LoD, int(map.sizeX - x - 1));
                int LoDY = std::min(LoD, int(map.sizeY - y - 1));
                int LoDZ = std::min(LoD, int(map.sizeZ - z - 1));
                float x_offset = x - (addedLeft ? 1.0 : 0.0);//-LoDIndex);
                float y_offset = y - (addedFront ? 1.0 : 0.0);//-LoDIndex);
                Vertex vertices[8] = {Vertex(x_offset       , y_offset       , z       , map(x       , y       , z       )),
                                      Vertex(x_offset + LoDX, y_offset       , z       , map(x + LoDX, y       , z       )),
                                      Vertex(x_offset + LoDX, y_offset + LoDY, z       , map(x + LoDX, y + LoDY, z       )),
                                      Vertex(x_offset       , y_offset + LoDY, z       , map(x       , y + LoDY, z       )),
                                      Vertex(x_offset       , y_offset       , z + LoDZ, map(x       , y       , z + LoDZ)),
                                      Vertex(x_offset + LoDX, y_offset       , z + LoDZ, map(x + LoDX, y       , z + LoDZ)),
                                      Vertex(x_offset + LoDX, y_offset + LoDY, z + LoDZ, map(x + LoDX, y + LoDY, z + LoDZ)),
                                      Vertex(x_offset       , y_offset + LoDY, z + LoDZ, map(x       , y + LoDY, z + LoDZ))
                                     };

                std::vector<Vector3> tempVerticesArray = this->computeMarchingCube(vertices, isolevel, useGlobalCoords, &colors, &normals);
                vertexArray.insert(vertexArray.end(), tempVerticesArray.begin(), tempVerticesArray.end());
            }
        }
    }
    if (outColors != nullptr)
        *outColors = colors;
    if (outNormals != nullptr)
        *outNormals = normals;
    return vertexArray;
}

void VoxelChunk::makeItFall()
{
    Matrix3<float> valuesAfterGravity(this->sizeX, this->sizeY, this->sizeZ);
    for (int x = 0; x < this->sizeX; x++) {
        for(int y = 0; y < this->sizeY; y++) {
            for(int z = 0; z < this->sizeZ - 1; z++) {
                if(this->voxelGroups(x, y, z) == 0 || this->voxelGroups(x, y, z+1) == 0)
                    continue;
                valuesAfterGravity.at(x, y, z) = -this->getVoxelValue(x, y, z) + this->getVoxelValue(x, y, z+1);
            }
            valuesAfterGravity.at(x, y, this->sizeZ - 1) = - this->getVoxelValue(x, y, this->sizeZ - 1) - 1.0;
        }
    }
    this->needRemeshing = true;
}
void VoxelChunk::letGravityMakeSandFall()
{
    Matrix3<float> transportMatrix(this->sizeX, this->sizeY, this->sizeZ);
    float sandLowerLimit = 0.0, sandUpperLimit = 1.0;
    for (int x = 0; x < this->sizeX; x++) {
        for(int y = 0; y < this->sizeY; y++) {
            int lastSandyUnsaturatedHeight = 0;
            float lastSandyUnsaturatedWetnessCoefficient = interpolation::smooth(lastSandyUnsaturatedHeight/(float)this->sizeZ);
            float lastSandySandUpperLimit = sandUpperLimit * remap(lastSandyUnsaturatedWetnessCoefficient, 0.f, 1.f, 0.2f, 1.f);
            for(int z = 0; z < this->sizeZ; z++) {
                float currentSandWetnessCoefficient = interpolation::smooth(z/(float)this->sizeZ);
                float currentSandUpperLimit = sandUpperLimit * remap(currentSandWetnessCoefficient, 0.f, 1.f, 0.2f, 1.f);
                float cellValue = this->parent->getVoxelValue(this->x + x, this->y + y, z) + transportMatrix.at(x, y, z);
                if(cellValue < sandLowerLimit) continue;
                else if (cellValue > currentSandUpperLimit) {
                    lastSandyUnsaturatedHeight = z + 1;
                    lastSandyUnsaturatedWetnessCoefficient = interpolation::smooth(lastSandyUnsaturatedHeight/(float)this->sizeZ);
                    lastSandySandUpperLimit = sandUpperLimit * remap(lastSandyUnsaturatedWetnessCoefficient, 0.f, 1.f, 0.2f, 1.f);
                    continue;
                } else //if (cellValue > currentSandUpperLimit)// If it's not dirt nor air
                {
                    // Make sand fall until the cell is empty or there's nowhere else to drop sand
                    while (cellValue > 0 && lastSandyUnsaturatedHeight != z) {
                        float lastSandyCellValue = this->parent->getVoxelValue(this->x + x, this->y + y, lastSandyUnsaturatedHeight) + transportMatrix.at(x, y, lastSandyUnsaturatedHeight);
                        float d_transported = std::min(cellValue, currentSandUpperLimit - lastSandyCellValue);
                        lastSandyCellValue += d_transported; // * 1.0001f;
                        cellValue -= d_transported;
                        transportMatrix.at(x, y, z) -= d_transported;
                        transportMatrix.at(x, y, lastSandyUnsaturatedHeight) += d_transported; // * 1.0001f;
                        if(lastSandyCellValue >= lastSandySandUpperLimit) {
                            lastSandyUnsaturatedHeight ++;
                            lastSandyUnsaturatedWetnessCoefficient = interpolation::smooth(lastSandyUnsaturatedHeight/(float)this->sizeZ);
                            lastSandySandUpperLimit = sandUpperLimit * remap(lastSandyUnsaturatedWetnessCoefficient, 0.f, 1.f, 0.2f, 1.f);
                        }
                    }
                }
            }
        }
    }
    this->applyModification(transportMatrix);
    this->needRemeshing = true;
}
void VoxelChunk::computeGroups()
{
    int currentMarker = 0;
    std::set<int> neighbors_coords;
    Matrix3<int> connected(this->sizeX, this->sizeY, this->sizeZ);
    Matrix3<int> labels(this->sizeX, this->sizeY, this->sizeZ, -1);
    Matrix3<float> voxelValues = this->getVoxelValues();
    for (size_t i = 0; i < voxelValues.size(); i++) {
        if (voxelValues[i] > 0.f) connected[i] = 1;
    }
    connected.raiseErrorOnBadCoord = false; // Allow to search on out-of-bounds
    connected.defaultValueOnBadCoord = 0; // But mark it as background

    for (size_t i = 0; i < connected.size(); i++) {
        if (connected.at(i) == 1) {
            currentMarker++;
            neighbors_coords.insert(i);
            while (!neighbors_coords.empty()) {
                int nx, ny, nz;
                size_t i_neighbor = *neighbors_coords.begin();
                neighbors_coords.erase(neighbors_coords.begin());
                labels.at(i_neighbor) = currentMarker;
                std::tie(nx, ny, nz) = connected.getCoord(i_neighbor);
                if (connected.at(i_neighbor) == 1) {
                    connected.at(i_neighbor) = 0;
                    for (int dx = -1; dx <= 1; dx++) {
                        for (int dy = -1; dy <= 1; dy++) {
                            for (int dz = -1; dz <= 1; dz++) {
                                if (connected(nx + dx, ny + dy, nz + dz) == 1) {
                                    neighbors_coords.insert(connected.getIndex(nx + dx, ny + dy, nz + dz));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    this->voxelGroups = labels;
}

void VoxelChunk::applyModification(Matrix3<float> modifications, Vector3 anchor)
{
    modifications.raiseErrorOnBadCoord = false;
    // If there is a modification and we did some "undo's" before,
    // erase the future matrices
    if (currentHistoryIndex < this->voxelsValuesStack.size()) {
        this->voxelsValuesStack.erase(this->voxelsValuesStack.begin() + currentHistoryIndex, this->voxelsValuesStack.end());
        this->voxelsValuesAnchorStack.erase(this->voxelsValuesAnchorStack.begin() + currentHistoryIndex, this->voxelsValuesAnchorStack.end());
    }
    this->currentHistoryIndex++;
    this->voxelsValuesStack.push_back(modifications);
    this->voxelsValuesAnchorStack.push_back(anchor);
    this->needRemeshing = true;
}

void VoxelChunk::undo()
{
    if (this->voxelsValuesStack.size() > 1) {
//        this->voxelsValuesStack.pop_back();
        this->currentHistoryIndex --;
        this->needRemeshing = true;
    }
}

void VoxelChunk::redo()
{
    if (this->currentHistoryIndex < this->voxelsValuesStack.size()) {
//        this->voxelsValuesStack.pop_back();
        this->currentHistoryIndex ++;
        this->needRemeshing = true;
    }
}


void VoxelChunk::display()
{
    this->mesh.display();
}

bool VoxelChunk::contains(Vector3 v) {
    return this->contains(v.x, v.y, v.z);
}

bool VoxelChunk::contains(float x, float y, float z) {
    return (this->x <= x && x < this->x + this->sizeX && this->y <= y && y < this->y + this->sizeY && 0 <= z && z < this->sizeZ);
}

float VoxelChunk::getVoxelValue(int x, int y, int z)
{
    float finalValue = 0.f;
    for (size_t i = 0; i < this->currentHistoryIndex; i++) {
        Vector3 anchor = this->voxelsValuesAnchorStack[i];
        finalValue += this->voxelsValuesStack[i].at(x - anchor.x, y - anchor.y, z - anchor.z);
    }
    return finalValue;
}

float VoxelChunk::getVoxelValue(Vector3 pos)
{
    return this->getVoxelValue(pos.x, pos.y, pos.z);
}

Matrix3<float> VoxelChunk::getVoxelValues()
{
    Matrix3<float> voxelValues(this->sizeX, this->sizeY, this->sizeZ);
    for (size_t i = 0; i < this->currentHistoryIndex; i++)
        voxelValues.add(this->voxelsValuesStack[i], this->voxelsValuesAnchorStack[i]);
    return voxelValues;
}
