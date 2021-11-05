#include "UnderwaterErosion.h"
#include "RockErosion.h"

UnderwaterErosion::UnderwaterErosion()
{

}
UnderwaterErosion::UnderwaterErosion(VoxelGrid* grid, int maxRockSize, float maxRockStrength, int rockAmount)
    : grid(grid), maxRockSize(maxRockSize), rockAmount(rockAmount), maxRockStrength(maxRockStrength)
{

}

void UnderwaterErosion::Apply()
{
    float starting_distance = pow(std::max(grid->sizeX, std::max(grid->sizeY, grid->sizeZ)) / 2.0, 2);
    starting_distance = sqrt(3 * starting_distance); // same as sqrt(x+y+z)
    int max_iter = 1000;
    int update_only_every = this->rockAmount / 10.0;
    for (int i = 0; i < this->rockAmount && max_iter > 0; i++)
    {
        int steps = 1000;
        Vector3 pos((.5 + .5 * (float)(random() / (float)RAND_MAX)) * starting_distance, (.5 + .5 * (float)(random() / (float)RAND_MAX)) * starting_distance, (.5 + .5 * (float)(random() / (float)RAND_MAX)) * starting_distance);
        Vector3 dir((float)(random() / (float)RAND_MAX) - .5, (float)(random() / (float)RAND_MAX) - .5, (float)(random() / (float)RAND_MAX) - .5);
        dir.normalize();
        RockErosion rock((float)(random() / (float)RAND_MAX) * this->maxRockSize, (float)(random() / (float)RAND_MAX) * this->maxRockStrength);

        bool touched = false;
        while (!touched && steps > 0) {
            steps --;
            pos += dir;
            Voxel* v = this->grid->getVoxel(pos.x, pos.y, pos.z);
            if (v != nullptr && v->type != TerrainTypes::AIR) {
                rock.Apply(v, false, (i % update_only_every == 0 || i == this->rockAmount));
                touched = true;
                std::cout << "Touche" << std::endl;
            }
            if (pos.norm() > 2 * starting_distance) {
                i--;
                max_iter --;
                break;
            }
        }
    }
}
