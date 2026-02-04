#ifndef CORALISLANDGENERATOR_H
#define CORALISLANDGENERATOR_H

#include "TerrainGen/Heightmap.h"
// #include "EnvObject/EnvPoint.h"
// #include "EnvObject/EnvCurve.h"
// #include "EnvObject/EnvArea.h"
#include "EnvObject/EnvironmentalScene.h"

class CoralIslandGenerator
{
public:
    CoralIslandGenerator();

    static GridF generate(GridF heights, float subsidence, float waterLevel, float coralMin, float maxCoralHeight, float verticalScale, float horizontalScale, float alpha);
    static std::vector<EnvObject*> envObjsFromFeatureMap(const GridV3& img, const Vector3 &terrainDimensions, std::shared_ptr<EnvironmentalScene> scene);

    // std::shared_ptr<EnvironmentalScene> scene;
};

#endif // CORALISLANDGENERATOR_H
