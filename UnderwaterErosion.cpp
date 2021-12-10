#include "Globals.h"
#include "UnderwaterErosion.h"
#include "RockErosion.h"
#include "BSpline.h"


UnderwaterErosion::UnderwaterErosion()
{

}
UnderwaterErosion::UnderwaterErosion(std::shared_ptr<VoxelGrid> grid, int maxRockSize, float maxRockStrength, int rockAmount)
    : grid(grid), maxRockSize(maxRockSize), rockAmount(rockAmount), maxRockStrength(maxRockStrength)
{

}

std::tuple<std::vector<std::vector<Vector3>>, std::vector<std::vector<Vector3>>>
UnderwaterErosion::Apply(std::shared_ptr<Vector3> startingPoint, std::shared_ptr<Vector3> originalDirection, int avoidMatter, float flowfieldFactor, float randomnessFactor, bool returnEvenLostRocks)
{
    /*for (std::shared_ptr<VoxelChunk>& vc : this->grid->chunks)
        vc->computeFlowfield();*/
    std::vector<std::vector<Vector3>> debugFinishingLines;
    std::vector<std::vector<Vector3>> debugFailingLines;
    float starting_distance = pow(std::max(grid->sizeX, std::max(grid->sizeY, grid->sizeZ))/2.0, 2);
    starting_distance = sqrt(3 * starting_distance); // same as sqrt(x+y+z)
    starting_distance *= 2.0; // Leave a little bit of gap
    int max_iter = 1000;
    int total_iterations = 0;
    int cpt = 0;
    for (int i = 0; i < this->rockAmount && max_iter > 0; i++)
    {
        float weaknessAgainstFlowfield = 1.0;
        cpt ++;
        int steps = 10 * starting_distance; // An estimation of how many step we need
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
            dir = (*originalDirection - pos).normalize();
        else if (startingPoint != nullptr && originalDirection != nullptr)
            dir = originalDirection->normalize();


        pos += Vector3(this->grid->getSizeX(), this->grid->getSizeY(), 0.0)/2.0;
        RockErosion rock(random_gen::generate(0.0, this->maxRockSize), random_gen::generate(0.0, this->maxRockStrength));
        std::vector<Vector3> coords;

        bool touched = false;
        while (!touched) {
            total_iterations ++;
            if (avoidMatter > 0) {
                for(int dist = 0; dist < avoidMatter; dist++) {
                    if (this->grid->contains(pos + dir * dist)) {
//                        dir += Vector3::random() * 0.05;
                        Vector3 flowfield = grid->getFlowfield(pos + dir * dist);
                        dir += flowfield * flowfieldFactor * (dist == 0 ? weaknessAgainstFlowfield : 1.0);
                        dir.normalize();
//                        if (dist == 0 && flowfield.norm() > 0.1)
//                            weaknessAgainstFlowfield *= .999f;//std::max(.95f, (1 - std::min(1.0f, flowfield.norm())));
                        break;
                    }
                }
            }
//            if (!matterIsClose) {
                dir += Vector3::random() * randomnessFactor;
                dir.normalize();
//            }
            steps --;
            coords.push_back(pos - Vector3(this->grid->getSizeX(), this->grid->getSizeY(), 0.0)/2.0);
            pos += dir;
            coords.push_back(pos - Vector3(this->grid->getSizeX(), this->grid->getSizeY(), 0.0)/2.0);
            std::shared_ptr<Voxel> v = this->grid->getVoxel(pos.x, pos.y, pos.z);
            if (v != nullptr && *v) {
                rock.Apply(v, false, false);
                touched = true;
                max_iter = 1000;
                debugFinishingLines.push_back(coords);
            }
            if (pos.norm() > 4 * starting_distance || steps <= 0 || pos.z < -100) {
                i--;
                max_iter --;
                if (returnEvenLostRocks)
                    debugFailingLines.push_back(coords);
                break;
            }
        }
    }
    std::cout << total_iterations << " iterations (" << cpt << " rocks lauched, " << debugFinishingLines.size() << " who hit)";
    std::cout << " check : " << debugFinishingLines.size() << "+" << debugFailingLines.size() << std::endl;
    grid->remeshAll();
    return std::make_tuple(debugFinishingLines, debugFailingLines);
}

std::vector<Vector3> UnderwaterErosion::CreateTunnel(std::shared_ptr<Vector3> startingPoint, std::shared_ptr<Vector3> endingPoint, int numberPoints, bool addingMatter)
{
    BSpline curve = BSpline(numberPoints); // Random curve
    curve.points[0].z = -1.0;
    curve.points[1].z =  0.9;
    curve.points[2].z = -1.0;
    BSpline width = BSpline(std::vector<Vector3>({
                                                     Vector3(0.0, 1.0),
                                                     Vector3(0.5, 0.5),
                                                     Vector3(1.0, 1.0)
                                                 }));

    float resolution = 0.01;
    std::vector<Vector3> coords;
    for (float i = 0; i < 1.0; i += resolution)
    {
        Vector3 pos = ((curve.getPoint(i) + 1.0) / 2.0) * Vector3(grid->sizeX, grid->sizeY, grid->sizeZ);
        coords.push_back(pos - Vector3(grid->sizeX/2.0, grid->sizeY/2.0, 0.0));
        float rockSize = width.getPoint(i).y * this->maxRockSize;
        RockErosion rock(random_gen::generate(0.0, rockSize), random_gen::generate(0.0, this->maxRockStrength));
        rock.Apply(grid->getVoxel(pos), addingMatter, false);
        coords.push_back((((curve.getPoint(i + resolution) + 1.0) / 2.0) * Vector3(grid->sizeX, grid->sizeY, grid->sizeZ)) - Vector3(grid->sizeX/2.0, grid->sizeY/2.0, 0.0));
    }
    grid->remeshAll();
    return coords;
}
