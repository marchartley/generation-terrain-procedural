#include "LayerBasedGrid.h"

#include "UnderwaterErosion.h"

#include "Shader.h"

#include "FastNoiseLit.h"
#include "Globals.h"

LayerBasedGrid::LayerBasedGrid() : LayerBasedGrid(10, 10, 10)
{

}
LayerBasedGrid::LayerBasedGrid(int nx, int ny, float nz)
    : sizeX(nx), sizeY(ny), sizeZ(nz)
{
    // Create and configure FastNoise object
    FastNoiseLite noise;
    noise.SetFrequency(0.01);
    noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
    noise.SetFractalType(FastNoiseLite::FractalType_FBm);
    noise.SetFractalLacunarity(2.0);
    noise.SetFractalGain(0.7);
    noise.SetFractalWeightedStrength(0.5);
    noise.SetFractalOctaves(10);

    std::vector<std::vector<std::vector<std::tuple<TerrainTypes, float>>>> data(nx);
    std::vector<std::vector<std::vector<std::shared_ptr<Voxel>>>> voxs(nx);
    for (int x = 0; x < nx; x++) {
        data[x] = std::vector<std::vector<std::tuple<TerrainTypes, float>>>(ny);
        voxs[x] = std::vector<std::vector<std::shared_ptr<Voxel>>>(ny);
        for (int y = 0; y < ny; y++) {
            data[x][y] = std::vector<std::tuple<TerrainTypes, float>>();
            voxs[x][y] = std::vector<std::shared_ptr<Voxel>>(int(nz));
            float remainingZ = nz;
            while (remainingZ > 0) {
                TerrainTypes type = static_cast<TerrainTypes>(int(random_gen::generate(0, TerrainTypes::LAST)));
                float height = 1.0; // noise.GetNoise((float)x, (float)y, remainingZ) + 1.0;
                remainingZ -= height;
                if (remainingZ < 0)
                    height += remainingZ; // Return to the exact max Z
                if (data[x][y].size() > 0 && std::get<0>(data[x][y][data[x][y].size() - 1]) == type) {
                    data[x][y][data[x][y].size() - 1] = std::make_pair<TerrainTypes&, float>(type, std::get<1>(data[x][y][data[x][y].size() - 1]) + height);
                }
                else {
                    data[x][y].push_back(std::make_tuple<TerrainTypes&, float&>(type, height));
                }
            }
            float current_height = 0.0;
            int current_level = 0;
            for (size_t i = 0; i < data[x][y].size(); i++) {
                for (int j = 0; j < std::get<1>(data[x][y][i]); j++) {
                    float iso = Voxel::isMatter[std::get<0>(data[x][y][i])] ? 1.0 : -1.0;
                    std::shared_ptr<Voxel> v(new Voxel(x, y, current_height, iso > 0 ? DIRT : AIR, 1.0, iso * std::get<1>(data[x][y][i]), nullptr));
                    voxs[x][y][current_level++] = v;
                    current_height += 1.0;
                }
            }
        }
    }
    this->layers = data;
    this->voxels = voxs;
}

void LayerBasedGrid::createMesh() {

//    this->computeGroups();

    std::vector<Vector3> voxelVertices;
    std::vector<Vector3> colors;
    this->applyToVoxels([&](std::shared_ptr<Voxel> v) -> void {
        if ((bool)*v) {
            // Add the vertices to the global mesh
            std::vector<Vector3> vertice = v->getMeshVertices();
            voxelVertices.insert(voxelVertices.end(), vertice.begin(), vertice.end());
            // Add the colors to each vertex
            int X = 6; // Start with 6 faces

            for(auto& n : v->neighbors)
                if (n.second && (bool)*n.second)
                    X--;    // Remove a face per neighbor
            X *= 6; // Multiply the number of face by the 6 vertex that defines it (2 triangles)
            for (int x = 0; x < X; x++) {
                colors.push_back(Vector3(random_gen::generate(), random_gen::generate(), random_gen::generate())); //(v->isOnGround ? 1.0 : 0.0), (v->isOnGround ? 0.0 : 1.0), 1.0));
//                        colors->push_back(HSVtoRGB((voxels[i][j][k]->group/((float)Voxel::voxelGroups.size()+1)), 1.0, 1.0));
            }
        }
    });

    this->mesh.colorsArray = colors;
    this->mesh.fromArray(voxelVertices);
    this->mesh.update();
}

void LayerBasedGrid::display() {
    this->mesh.display();
}
    /*
    // Create and configure FastNoise object
    FastNoiseLite noise;
    noise.SetFrequency(0.01);
    noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
    noise.SetFractalType(FastNoiseLite::FractalType_FBm);
    noise.SetFractalLacunarity(2.0);
    noise.SetFractalGain(0.7);
    noise.SetFractalWeightedStrength(0.5);
    noise.SetFractalOctaves(10);

    std::vector<std::vector<std::vector<std::tuple<TerrainTypes, float>>>> data(nx);
    std::vector<std::vector<std::vector<std::shared_ptr<Voxel>>>> voxs(nx);
    for (int x = 0; x < nx; x++) {
        data[x] = std::vector<std::vector<std::tuple<TerrainTypes, float>>>(ny);
        voxs[x] = std::vector<std::vector<std::shared_ptr<Voxel>>>(ny);
        for (int y = 0; y < ny; y++) {
            data[x][y] = std::vector<std::tuple<TerrainTypes, float>>();
            voxs[x][y] = std::vector<std::shared_ptr<Voxel>>(int(nz));
            float remainingZ = nz;
            while (remainingZ > 0) {
                TerrainTypes type = static_cast<TerrainTypes>(int(random_gen::generate(0, TerrainTypes::LAST)));
                float height = noise.GetNoise((float)x, (float)y, remainingZ) + 1.0;
                remainingZ -= height;
                if (remainingZ < 0)
                    height += remainingZ; // Return to the exact max Z
                if (data[x][y].size() > 0 && std::get<0>(data[x][y][data[x][y].size() - 1]) == type) {
                    data[x][y][data[x][y].size() - 1] = std::make_pair<TerrainTypes&, float>(type, std::get<1>(data[x][y][data[x][y].size() - 1]) + height);
                }
                else {
                    data[x][y].push_back(std::make_tuple<TerrainTypes&, float&>(type, height));
                }
            }
            for (int i = 0; i < nz; i++) {
                std::shared_ptr<Voxel> v(new Voxel(x, y, i, DIRT, 10.0, 1.0, nullptr));
                voxs[x][y][i] = v;
            }
        }
    }
    this->layers = data;
    this->voxels = voxs;
}

*/
