#include "Globals.h"
#include "UnderwaterErosion.h"
#include "RockErosion.h"



UnderwaterErosion::UnderwaterErosion()
{

}
UnderwaterErosion::UnderwaterErosion(VoxelGrid* grid, int maxRockSize, float maxRockStrength, int rockAmount)
    : grid(grid), maxRockSize(maxRockSize), rockAmount(rockAmount), maxRockStrength(maxRockStrength)
{

}

std::vector<std::vector<Vector3>> UnderwaterErosion::Apply(Vector3* startingPoint, Vector3* originalDirection, int avoidMatter)
{
    std::vector<std::vector<Vector3>> debugLines;
    float starting_distance = pow(std::max(grid->sizeX, std::max(grid->sizeY, grid->sizeZ)) / 2.0, 2);
    starting_distance = sqrt(3 * starting_distance); // same as sqrt(x+y+z)
    int max_iter = 1000;
    for (int i = 0; i < this->rockAmount && max_iter > 0; i++)
    {
        int steps = 1000;
        Vector3 pos;
        if (startingPoint == nullptr) {
            pos = Vector3(random_gen::generate(-1.0, 1.0), random_gen::generate(-1.0, 1.0), random_gen::generate(0.0, 1.0));
            pos.normalize();
            pos *= starting_distance;
        } else {
            pos = *startingPoint;
        }
        Vector3 dir = Vector3::random();
        if (startingPoint == nullptr && originalDirection != nullptr)
            dir = (originalDirection - pos).normalize();
        else if (startingPoint != nullptr && originalDirection != nullptr)
            dir = originalDirection->normalize();

        pos += Vector3(this->grid->getSizeX(), this->grid->getSizeY(), 0.0)/2.0;
        RockErosion rock(random_gen::generate(0.0, this->maxRockSize), random_gen::generate(0.0, this->maxRockStrength));
        std::vector<Vector3> coords;

        bool touched = false;
        while (!touched && steps > 0) {
            dir += Vector3::random() * 0.1;
            dir.normalize();
            // Try to change a little bit the direction if there is matter ahead
            if (avoidMatter > 0) {
                for(int tries = 0; tries < avoidMatter/5; tries ++) {
                    for(int dist = 1; dist < avoidMatter; dist++) {
                        Voxel* v = this->grid->getVoxel(pos + dir * dist);
                        if (v != nullptr && *v) {
                            dir += Vector3::random() * 0.3;
                            dir.normalize();
                            break;
                        }
                    }
                }
            }
            steps --;
            pos += dir;
            coords.push_back(pos - Vector3(this->grid->getSizeX(), this->grid->getSizeY(), 0.0)/2.0);
            Voxel* v = this->grid->getVoxel(pos.x, pos.y, pos.z);
            if (v != nullptr && *v) {
                rock.Apply(v, false, false);
                touched = true;
                max_iter = 1000;
                debugLines.push_back(coords);
            }
            if (pos.norm() > 2 * starting_distance) {
                i--;
                max_iter --;
                break;
            }
        }
    }
    grid->remeshAll();
    return debugLines;
}
