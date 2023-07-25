#include "SpheroidalWeathering.h"

SpheroidalWeathering::SpheroidalWeathering(std::shared_ptr<VoxelGrid> voxelGrid)
    : voxelGrid(voxelGrid)
{
    if (voxelGrid) {
        voxels = voxelGrid->getVoxelValues().binarize();
        decimation = Matrix3<float>(voxels.getDimensions(), 1.f);
    }
}

void SpheroidalWeathering::applyErosion()
{
    // "function f [...] can scale or clamp the curvature estimate to produce different effects"
    std::function<float(float)> curvatureTransform = [](float someValue) -> float {
        return someValue - .1f;
    };

    float bubbleRadius = 1.f;

    if (decimation.size() == 0) {
        voxels = voxelGrid->getVoxelValues().binarize();
        decimation = Matrix3<float>(voxels.getDimensions(), 1.f);
    }

    Matrix3<float> airPercents = Matrix3<float>(voxels.getDimensions());
    Matrix3<float> layerDurability = Matrix3<float>(voxels.getDimensions());
    FastNoiseLite noise;
    noise.SetFractalType(FastNoiseLite::FractalType_FBm);
    float minDurability = 0.f, maxDurability = 1.f;
    for (int x = 0; x < voxels.sizeX; x++) {
        for (int y = 0; y < voxels.sizeY; y++) {
            for (int z = 0; z < voxels.sizeZ; z++) {
                float noiseVal = noise.GetNoise((float) z * 5.f, x*y/100.f);
                noiseVal = (noiseVal + 1.f) * .5f;
                layerDurability.at(x, y, z) = (noiseVal * (maxDurability - minDurability)) + minDurability;
                if (z > voxels.sizeZ * .8) {
                    float dz = (voxels.sizeZ - z) / (voxels.sizeZ * .2);
                    layerDurability.at(x, y, z) *= interpolation::wyvill(1.f - dz);
                } else {
                    layerDurability.at(x, y, z) *= interpolation::wyvill(1.f - z / (voxels.sizeZ * .8));
                }
            }
        }
    }

    for (int _ = 0; _ < 1; _++) {

        for (auto& v : voxels)
            v = (v < .25f ? 0.f : (v < .75f ? .5f : 1.f));
    //    voxels = voxels.binarize();
        auto surface = voxels.binarize() - voxels.binarize().erode();

        for (int x = 0; x < surface.sizeX; x++) {
            for (int y = 0; y < surface.sizeY; y++) {
                for (int z = 0; z < 3; z++) {
                    surface.at(x, y, z) = 0;
                }
            }
        }
        std::cout << "Surface voxels: " << surface.sum() << std::endl;

        for (int x = 0; x < voxels.sizeX; x++) {
            for (int y = 0; y < voxels.sizeY; y++) {
                for (int z = 0; z < voxels.sizeZ; z++) {
                    Vector3 pos = Vector3(x, y, z);
                    if (surface.at(pos) == 0) continue;
                    int nbRocks = 0;
                    int nbAir = 0;
                    for (float dx = -bubbleRadius; dx <= bubbleRadius; dx++) {
                        for (float dy = - bubbleRadius; dy <= bubbleRadius; dy++) {
                            for (float dz = -bubbleRadius; dz <= bubbleRadius; dz++) {
                                Vector3 offset = Vector3(dx, dy, dz);
                                if (offset.norm2() < bubbleRadius * bubbleRadius) continue;
                                if (voxels.at(pos + offset) == 0)
                                    nbAir++;
                                else
                                    nbRocks++;
                            }
                        }
                    }
                    airPercents.at(pos) = float(nbAir) / float(nbRocks + nbAir);
                    decimation.at(pos) = decimation.at(pos) - (curvatureTransform(airPercents.at(pos)) / layerDurability.at(pos));
//                    decimation.at(pos) = decimation.at(pos) - (curvatureTransform(1.f - airPercents.at(pos)) / layerDurability.at(pos)); // Not sure the 2 calls are useful...

    //                std::cout << airPercents.at(pos) << " of air -> decimation = " << decimation.at(pos) << " (layerDurability is " << layerDurability.at(pos) << ")" << std::endl;
                    if (decimation.at(pos) <= 0) {
    //                    std::cout << "From " << voxels.at(pos);
                        _erode(voxels, pos);
    //                    std::cout << " to " << voxels.at(pos) << std::endl;
                    } else {

                    }
                }
            }
        }
        this->_deposition(voxels);
    }
    voxelGrid->setVoxelValues(voxels - .25f);
}

void SpheroidalWeathering::_erode(Matrix3<float> &voxels, const Vector3 &pos)
{
    if (voxels.at(pos) <= .75f) { // Is colluvium
//        voxels.at(pos) = 0.f;
    } else { // Is rock
        voxels.at(pos) = .5f;

        // More to do to create deposition
    }
}

void SpheroidalWeathering::_deposition(Matrix3<float> &voxels)
{
    float talusHeight = 3; // std::atan(deg2rad(30));

//    Matrix3<float> _voxels = voxels;

    for (int x = 0; x < voxels.sizeX; x++) {
        for (int y = 0; y < voxels.sizeY; y++) {
            for (int z = 0; z < voxels.sizeZ; z++) {
                if (voxels.at(x, y, z) < .25f || .75f < voxels.at(x, y, z)) continue;
                float slope = 1000;
                Vector3 newPosition(x, y, z);
                while(slope > talusHeight) {
                    Vector3 bestDirection;
                    for (int dx = -1; dx <= 1; dx++) {
                        for (int dy = -1; dy <= 1; dy++) {
                            int depth = 0;
                            while (voxels.at(newPosition + Vector3(dx, dy, depth)) < .25f && z + depth > 0) {
                                depth--;
                            }
                            if (depth < bestDirection.z) {
                                bestDirection = Vector3(dx, dy, depth + 1);
                            }
                        }
                    }
                    slope = -bestDirection.z;
                    newPosition += bestDirection;
                }
                voxels.at(x, y, z) = 0;
                voxels.at(newPosition) = .5f;
            }
        }
    }
}
