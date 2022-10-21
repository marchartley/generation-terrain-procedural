#ifndef UNDERWATEREROSION_H
#define UNDERWATEREROSION_H

#include "TerrainGen/VoxelGrid.h"

#include <tuple>
#include "Utils/BSpline.h"
#include "Karst/KarstHoleProfile.h"
#include "Karst/KarstHole.h"

class UnderwaterErosion
{
public:
    UnderwaterErosion();
    UnderwaterErosion(std::shared_ptr<VoxelGrid> grid, int maxRockSize, float maxRockStrength, int rockAmount);

//    std::vector<std::vector<Vector3>> Apply(int avoidMatter = -1);
    std::tuple<std::vector<std::vector<Vector3>>, std::vector<std::vector<Vector3>>> Apply(Vector3 startingPoint = Vector3(false),
                                                                                           Vector3 originalDirection = Vector3(false),
                                                                                           float randomnessFactor = 0.05,
                                                                                           bool fallFromSky = false,
                                                                                           float gravity = 1.f,
                                                                                           float bouncingCoefficient = .5f,
                                                                                           float bounciness = .5f,
                                                                                           float minSpeed = .1f,
                                                                                           float maxSpeed = 1.f,
                                                                                           float maxCapacityFactor = 1.f,
                                                                                           float erosion = 1.f,
                                                                                           float deposit = 1.f,
                                                                                           float matterDensity = 1000.f,
                                                                                           float materialImpact = 0.f,
                                                                                           float airFlowfieldRotation = 0.f,
                                                                                           float waterFlowfieldRotation = 0.f,
                                                                                           float airForce = 1.f,
                                                                                           float waterForce = 1.f);

    std::vector<Vector3> CreateTunnel(int numberPoints = 2, bool addingMatter = false, bool applyChanges = true, KarstHolePredefinedShapes startingShape = SOLUBLE_BED, KarstHolePredefinedShapes endingShape = KEYHOLE);
    std::vector<Vector3> CreateTunnel(BSpline path, bool addingMatter = false, bool usingSpheres = true, bool applyChanges = true, KarstHolePredefinedShapes startingShape = SOLUBLE_BED, KarstHolePredefinedShapes endingShape = KEYHOLE);
    std::vector<std::vector<Vector3>> CreateMultipleTunnels(std::vector<BSpline> paths, bool addingMatter = false, bool usingSpheres = true, bool applyChanges = true, KarstHolePredefinedShapes startingShape = SOLUBLE_BED, KarstHolePredefinedShapes endingShape = KEYHOLE);
    std::vector<Vector3> CreateCrack(Vector3 start, Vector3 end, bool applyChanges = true);

    // With a predefined tunnel :
    std::vector<Vector3> CreateTunnel(KarstHole& tunnel, bool addingMatter = false, bool applyChanges = true);


    std::shared_ptr<VoxelGrid> grid;
    int maxRockSize;
    int rockAmount;
    float maxRockStrength;

};

#endif // UNDERWATEREROSION_H
