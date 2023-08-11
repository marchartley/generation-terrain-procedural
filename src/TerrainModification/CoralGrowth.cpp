#include "CoralGrowth.h"

CoralGrowth::CoralGrowth()
{
    volume = GridF(terrainSize);
    coralArea = GridF(terrainSize);
    highErosionArea = GridF(terrainSize);
    highDepositArea = GridF(terrainSize);

    for (int x = 0; x < terrainSize.x; x++) {
        for (int y = 0; y < terrainSize.y; y++) {
            for (int z = 0; z < terrainSize.z; z++) {
                if (5 < y && y < 10 && z < 5) {
                    coralArea.at(x, y, z) = 1.f;
                } else {
                    coralArea.at(x, y, z) = -.01f;
                }
            }
        }
    }
    volume = coralArea;
}

void CoralGrowth::step()
{
    highErosionArea.reset();
    highDepositArea.reset();

    GridF surface = (1.f - volume.binarize()).toDistanceMap();
    GridV3 normals = surface.gradient();
    for (auto& n : normals)
        n.normalize();
    for (auto& val : surface)
        val = (val == 0.f ? 1.f : 0.f);

    // Find most exposed and least exposed areas
    for (int x = 0; x < terrainSize.x; x++) {
        for (int y = 0; y < terrainSize.y; y++) {
            for (int z = 0; z < terrainSize.z; z++) {
                Vector3 pos(x, y, z);
                if (!surface.at(pos))
                    continue;
                float shadeDistance = 0.f;
                bool shaded = false;
                Vector3 checkingPos = pos;
                while (Vector3::isInBox(checkingPos, Vector3(), terrainSize) && !shaded) {
                    checkingPos = checkingPos - waterflow + Vector3::random(0.5f);
                    if (volume.at(checkingPos))
                        shaded = true;
                }

                if (shaded) {
                    shadeDistance = (checkingPos - pos).norm() / terrainSize.x;
                    highDepositArea.at(pos) = (1.f - shadeDistance);
                }  else {
                    highErosionArea.at(pos) = std::abs(1.f - normals.at(pos).dot(waterflow));
                }

                if (z > waterHeight - 10) {
                    highErosionArea.at(pos) *= 2.f;
                }
            }
        }
    }

    // Change the volume depending on erosion/deposition
    float totalErosion = highErosionArea.sum();
    float totalDeposit = highDepositArea.sum();

    float ratio = totalErosion / totalDeposit;

    highDepositArea *= ratio;

    GridF modifications = (highDepositArea - highErosionArea).meanSmooth(5, 5, 5);

    volume += modifications;
}
