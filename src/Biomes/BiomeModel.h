#ifndef BIOMEMODEL_H
#define BIOMEMODEL_H

class BiomeModel;

#include "Biomes/BiomeInstance.h"
#include "Utils/json.h"
#include "DataStructure/Vector3.h"
#include "Utils/ShapeCurve.h"

#include <cmath>

enum BIOME_SLOPE_TYPE {
    MOUNTAIN     = 0,
    CURVE        = 1,
    FLAT         = 2,
    LINEAR_SLOPE = 3,
    UNDEFINED    = 4
};

struct Probability {
    float mean;
    float variance;
    Probability(float mean = 0.f, float variance = 0.f) : mean(mean), variance(variance) {}
    float randomValue() {
        return random_gen::generate_normal(mean, std::sqrt(variance));
    }
    nlohmann::json toJson() {
        return nlohmann::json({{"mean", mean}, {"vari", variance}});
    }
    static Probability fromJson(nlohmann::json content) {
        Probability proba;
        if (content.contains("mean")) {
            proba.mean = content.at("mean").get<float>();
            proba.variance = content.at("vari").get<float>();
        } else {
            proba.mean = content.get<float>();
        }
        return proba;
    }
};

struct BiomeModelChild {
    std::shared_ptr<BiomeModel> model;

    Probability probaQuantity;
    Probability probaAppearence;
    Probability probaSize;
    int priorityOffset;
};
class BiomeModel : public std::enable_shared_from_this<BiomeModel>
{
public:
    BiomeModel();

    static BiomeModel fromJson(nlohmann::json json_content);

    nlohmann::json toJson();
    std::shared_ptr<BiomeInstance> createInstance(const Vector3& initialPosition = Vector3(500, 500, 0),
                                 ShapeCurve initialArea = ShapeCurve({Vector3(0, 0, 0),
                                                                     Vector3(0, 1000, 0),
                                                                     Vector3(1000, 1000, 0),
                                                                     Vector3(1000, 0, 0)}));

    std::shared_ptr<BiomeModel> clone();

//    BSpline depthShape;
    BIOME_SLOPE_TYPE slopeType = BIOME_SLOPE_TYPE::UNDEFINED;
    float minDepth, maxDepth;
    std::string textureClass;
    std::string modelName;
    int priortyOffset = 0;
    Probability idealSize;
    std::vector<BiomeModelChild> modelChildren;
};

#endif // BIOMEMODEL_H
