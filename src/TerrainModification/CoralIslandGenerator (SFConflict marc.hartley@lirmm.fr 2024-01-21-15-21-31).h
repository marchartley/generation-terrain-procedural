#ifndef CORALISLANDGENERATOR_H
#define CORALISLANDGENERATOR_H

#include "TerrainGen/Heightmap.h"
#include "EnvObject/EnvObject.h"

class CoralIslandGenerator
{
public:
    CoralIslandGenerator();

    static GridF generate(GridF heights, float subsidence, float waterLevel, float coralMin, float maxCoralHeight, float verticalScale, float horizontalScale, float alpha);
    static std::vector<EnvObject*> envObjsFromFeatureMap(const GridV3& img, const Vector3 &terrainDimensions);
};

#endif // CORALISLANDGENERATOR_H
