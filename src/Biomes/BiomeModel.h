#ifndef BIOMEMODEL_H
#define BIOMEMODEL_H

class BiomeModel;

#include "Biomes/BiomeInstance.h"
#include "Utils/json.h"
#include "DataStructure/Vector3.h"
#include "Utils/ShapeCurve.h"

class BiomeModel
{
public:
    BiomeModel();

    static BiomeModel fromJson(nlohmann::json json_content);
    BiomeInstance createInstance(Vector3 initialPosition = Vector3(500, 500, 0),
                                 ShapeCurve initialArea = ShapeCurve({Vector3(0, 0, 0),
                                                                     Vector3(0, 1000, 0),
                                                                     Vector3(1000, 1000, 0),
                                                                     Vector3(1000, 0, 0)}));
    BSpline depthShape;
    std::string textureClass;
    std::string modelName;
    std::vector<BiomeModel> modelChildren;
    std::vector<BiomeInstance> instanceChildren;

    // Add the parameters in here
};

#endif // BIOMEMODEL_H
