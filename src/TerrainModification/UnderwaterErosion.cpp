#include "Utils/Globals.h"
#include "TerrainModification/UnderwaterErosion.h"
#include "TerrainModification/RockErosion.h"
#include "Utils/BSpline.h"


UnderwaterErosion::UnderwaterErosion()
{

}
UnderwaterErosion::UnderwaterErosion(std::shared_ptr<VoxelGrid> grid, int maxRockSize, float maxRockStrength, int rockAmount)
    : grid(grid), maxRockSize(maxRockSize), rockAmount(rockAmount), maxRockStrength(maxRockStrength)
{

}

void retroChangeFlowfield(std::vector<Vector3>& coords, std::vector<Vector3>& dirs, std::shared_ptr<VoxelGrid> grid)
{
    if (coords.size() == 0)
        return;
    Vector3& impactZone = coords[coords.size() - 1];
    for (size_t i = 0; i < coords.size(); i++)
    {
        Vector3& coord = coords[i];
        Vector3& dir = dirs[i];
        float alpha_effect = 0.5 * ((coord - impactZone).norm2() < 20.0 ? -2.0 : 1.0); // Inverse and double if it's a choc (last coord)
        grid->affectFlowfieldAround(coord, dir * alpha_effect, 3);
//        grid->affectFlowfieldAround(coord, alpha_effect, 10);
    }
}

//void updateParticle(Vector3& pos, Vector3& dir)

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
    std::vector<std::tuple<Vector3, Vector3>> usefulStartingPositions;
    for (int i = 0; i < this->rockAmount && max_iter > 0; i++)
    {
        std::vector<Vector3> allCoords;
        std::vector<Vector3> allDirs;
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
        /*if (startingPoint == nullptr && originalDirection != nullptr)
            dir = (*originalDirection - pos).normalize();
        else */if (/*startingPoint != nullptr &&*/ originalDirection != nullptr)
            dir = originalDirection->normalize();
//        dir += Vector3::random() * .1;


//        pos += Vector3(this->grid->getSizeX(), this->grid->getSizeY(), 0.0)/2.0;
        RockErosion rock(random_gen::generate(0.0, this->maxRockSize), random_gen::generate(0.0, this->maxRockStrength));
        std::vector<Vector3> coords;

        Vector3 lastSavedPos = pos - 10.0;
        bool firstPosSave = true;
        bool touched = false;
        bool hasBeenAtLeastOnceInside = false;
        float dist = 0;
        std::tuple<Vector3, Vector3> possibleUsefulPosition = std::make_tuple(pos, dir);
        while (!touched) {
            total_iterations ++;
            if (this->grid->contains(pos + dir * dist)) {
                Vector3 flowfield = grid->getFlowfield(pos + dir * dist);
                float dirDotGrad = 1.0; // dir.dot(flowfield);
                dir += flowfield * dirDotGrad * flowfieldFactor * (dist == 0 ? weaknessAgainstFlowfield : 1.0);
                dir.normalize(); // Maybe to remove, who knows...
                hasBeenAtLeastOnceInside = true;
//                break;
            }
            else {
                dir += Vector3::random() * randomnessFactor;
                dir.normalize();
            }
            steps --;
            pos += dir.normalized();
            if (firstPosSave) {
                firstPosSave = false;
                coords.push_back((pos - dir.normalized())); // - Vector3(this->grid->getSizeX(), this->grid->getSizeY(), 0.0)/2.0);
                lastSavedPos = pos;
            }
            else if ((lastSavedPos - pos).norm2() > 1.0) {
                coords.push_back(pos); // - Vector3(this->grid->getSizeX(), this->grid->getSizeY(), 0.0)/2.0);
                coords.push_back(pos); // - Vector3(this->grid->getSizeX(), this->grid->getSizeY(), 0.0)/2.0);
                lastSavedPos = pos;
            }
            allCoords.push_back(pos);
            allDirs.push_back(dir);
            if (this->grid->getVoxelValue(pos) > 0.0) { // Hit a wall
                if (hasBeenAtLeastOnceInside) {
                    rock.Apply(this->grid, pos, false, false);
                    coords.push_back(pos); // - Vector3(this->grid->getSizeX(), this->grid->getSizeY(), 0.0)/2.0);
                    debugFinishingLines.push_back(coords);
                    // Set a little bit the flowfield to avoid other particles going in same wall (increase of pressure)
//                    retroChangeFlowfield(allCoords, allDirs, this->grid);
//                    this->grid->affectFlowfieldAround(pos, dir.normalized() * -1.0f, 3);
                    usefulStartingPositions.push_back(possibleUsefulPosition);
                } else {
//                    i --;
                }
                touched = true;
                max_iter = 1000;
            }
            else if ((!hasBeenAtLeastOnceInside && (pos + dir).norm2() > pos.norm2()) || pos.norm2() > 16 * starting_distance * starting_distance || steps <= 0 || pos.z < -100 || pos.z > this->grid->getSizeZ()) {
                // Failed to hit, is now gone
                i--;
                max_iter --;
                if (returnEvenLostRocks && hasBeenAtLeastOnceInside) {
                    coords.push_back(pos); // - Vector3(this->grid->getSizeX(), this->grid->getSizeY(), 0.0)/2.0);
                    debugFailingLines.push_back(coords);
                }
                break;
            } else { // Particle still running
                // Set a little bit the flowfield to bring other particles in the same path
//                this->grid->affectFlowfieldAround(pos, dir.normalized() * 1.0f, 3);
            }
        }
    }
    std::cout << total_iterations << " iterations (" << cpt << " rocks lauched, " << debugFinishingLines.size() << " who hit)" << std::flush;
    grid->remeshAll();
    std::cout << " check : " << debugFinishingLines.size() << "+" << debugFailingLines.size() << std::endl;
    return std::make_tuple(debugFinishingLines, debugFailingLines);
}


std::vector<Vector3> UnderwaterErosion::CreateTunnel(int numberPoints, bool addingMatter)
{
    BSpline curve = BSpline(numberPoints); // Random curve
    for (Vector3& coord : curve.points)
        coord = ((coord + 1.0) / 2.0) * Vector3(grid->sizeX, grid->sizeY, grid->sizeZ);
    return CreateTunnel(curve, addingMatter);
}
std::vector<Vector3> UnderwaterErosion::CreateTunnel(BSpline path, bool addingMatter)
{
    Matrix3<float> erosionMatrix(this->grid->sizeX, this->grid->sizeY, this->grid->sizeZ);
    BSpline width = BSpline(std::vector<Vector3>({
                                                     Vector3(0.0, 1.0),
                                                     Vector3(0.5, 0.5),
                                                     Vector3(1.0, 1.0)
                                                 }));

    float resolution = 1.0 / path.length();
    std::vector<Vector3> coords;
    for (float i = 0; i < 1.0; i += resolution)
    {
        Vector3 pos = (path.getPoint(i));
        coords.push_back(pos); // - Vector3(grid->sizeX/2.0, grid->sizeY/2.0, 0.0));
        float rockSize = width.getPoint(i).y * this->maxRockSize;
        RockErosion rock(random_gen::generate(0.0, rockSize), random_gen::generate(0.0, this->maxRockStrength));
        erosionMatrix = rock.computeErosionMatrix(erosionMatrix, pos, addingMatter);
        coords.push_back(path.getPoint(i + resolution)); // - Vector3(grid->sizeX/2.0, grid->sizeY/2.0, 0.0));
    }
    grid->applyModification(erosionMatrix);
    grid->remeshAll();
    return coords;
}

std::vector<std::vector<Vector3> > UnderwaterErosion::CreateMultipleTunnels(std::vector<BSpline> paths, bool addingMatter)
{
    Matrix3<float> erosionMatrix(this->grid->sizeX, this->grid->sizeY, this->grid->sizeZ);
    BSpline width = BSpline(std::vector<Vector3>({
                                                     Vector3(0.0, 1.0),
                                                     Vector3(0.5, 0.5),
                                                     Vector3(1.0, 1.0)
                                                 }));

//    float resolution = 0.01;
    std::vector<std::vector<Vector3>> allCoords;
    for (BSpline& path : paths) {
        float resolution = 0.10 / path.length();
        std::vector<Vector3> coords;
        for (float i = 0; i < 1.0; i += resolution)
        {
            Vector3 pos = path.getPoint(i);
            coords.push_back(pos); // - Vector3(grid->sizeX/2.0, grid->sizeY/2.0, 0.0));
            float rockSize = width.getPoint(i).y * this->maxRockSize;
            RockErosion rock(random_gen::generate(0.0, rockSize), random_gen::generate(0.0, this->maxRockStrength));
            erosionMatrix = rock.computeErosionMatrix(erosionMatrix, pos, addingMatter);
            coords.push_back(path.getPoint(i + resolution)); // - Vector3(grid->sizeX/2.0, grid->sizeY/2.0, 0.0));
        }
        allCoords.push_back(coords);
    }
    grid->applyModification(erosionMatrix);
    grid->remeshAll();
    return allCoords;
}

