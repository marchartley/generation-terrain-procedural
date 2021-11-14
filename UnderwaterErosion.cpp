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

std::vector<std::vector<Vector3>> UnderwaterErosion::Apply()
{
    std::vector<std::vector<Vector3>> debugLines;
    float starting_distance = pow(std::max(grid->sizeX, std::max(grid->sizeY, grid->sizeZ)) / 2.0, 2);
    starting_distance = sqrt(3 * starting_distance); // same as sqrt(x+y+z)
    int max_iter = 1000;
//    int update_only_every = this->rockAmount / 10.0;
    for (int i = 0; i < this->rockAmount && max_iter > 0; i++)
    {
        int steps = 1000;
        Vector3 pos(random_gen::generate(-1.0, 1.0), random_gen::generate(-1.0, 1.0), random_gen::generate(0.0, 1.0));

        pos.normalize();
        pos *= starting_distance;
        pos += Vector3(this->grid->getSizeX(), this->grid->getSizeY(), this->grid->getSizeZ())/2.0;
        Vector3 dir = Vector3::random();
        RockErosion rock(random_gen::generate(0.0, this->maxRockSize), random_gen::generate(0.0, this->maxRockStrength));
        std::vector<Vector3> coords;

        bool touched = false;
        while (!touched && steps > 0) {
            dir += Vector3::random() * 0.1;
            dir.normalize();
            steps --;
            pos += dir;
            coords.push_back(pos);
            Voxel* v = this->grid->getVoxel(pos.x, pos.y, pos.z);
            if (v != nullptr && v->type != TerrainTypes::AIR) {
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
std::vector<std::vector<Vector3>> UnderwaterErosion::Apply(Vector3 startingPoint)
{
    std::vector<std::vector<Vector3>> debugLines;
    float starting_distance = pow(std::max(grid->sizeX, std::max(grid->sizeY, grid->sizeZ)) / 2.0, 2);
    starting_distance = sqrt(3 * starting_distance); // same as sqrt(x+y+z)
    int max_iter = 1000;
//    int update_only_every = this->rockAmount / 10.0;
    for (int i = 0; i < this->rockAmount && max_iter > 0; i++)
    {
        int steps = 1000;
//        Vector3 pos(random_gen::generate(-1.0, 1.0), random_gen::generate(-1.0, 1.0), random_gen::generate(0.0, 1.0));
        Vector3 pos = startingPoint;
//        pos.normalize();
//        pos *= starting_distance;
        pos += Vector3(this->grid->getSizeX(), this->grid->getSizeY(), this->grid->getSizeZ())/2.0;
        Vector3 dir = Vector3::random();
        RockErosion rock(random_gen::generate(0.0, this->maxRockSize), random_gen::generate(0.0, this->maxRockStrength));
        std::vector<Vector3> coords;

        bool touched = false;
        while (!touched && steps > 0) {
            dir += Vector3::random() * 0.1;
            dir.normalize();
            steps --;
            pos += dir;
            coords.push_back(pos);
            Voxel* v = this->grid->getVoxel(pos.x, pos.y, pos.z);
            if (v != nullptr && v->type != TerrainTypes::AIR) {
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
