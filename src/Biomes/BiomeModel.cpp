#include "BiomeModel.h"
#include "Utils/ConstraintsSolver.h"
#include "Utils/Voronoi.h"

BiomeModel::BiomeModel()
{

}

BiomeModel BiomeModel::fromJson(nlohmann::json json_content)
{
    BiomeModel model;
    if (json_content.contains("model-name") || json_content.contains("class"))
    {
        if (json_content.contains("model-name"))
            model.modelName = json_content.at("model-name").get<std::string>();
        else if (json_content.contains("class"))
            model.modelName = json_content.at("class").get<std::string>();

        // Extract parameters
        if (json_content.contains("params")) {
            auto params = json_content.at("params");
            if (params.contains("depth")) {
                auto depthInfo = params.at("depth");
                float minDepth = depthInfo.at("min").get<float>();
                float maxDepth = depthInfo.at("max").get<float>();
                if (depthInfo.contains("slope")) {
                    std::string slopeType = depthInfo.at("slope").get<std::string>();
                    if (slopeType == "mountain") {
                        model.depthShape = BSpline({
                            Vector3(0.0, maxDepth),
                            Vector3(0.5, minDepth),
                            Vector3(1.0, maxDepth)
                        });
                    } else if (slopeType == "curve") {
                        model.depthShape = BSpline({
                            Vector3(0.0, minDepth),
                            Vector3(0.5, maxDepth),
                            Vector3(1.0, minDepth)
                        });
                    } else if (slopeType == "flat") {
                        model.depthShape = BSpline({
                            Vector3(0.0, minDepth),
                            Vector3(1.0, minDepth)
                        });
                    } else if (slopeType == "linear") {
                        model.depthShape = BSpline({
                            Vector3(0.0, minDepth),
                            Vector3(1.0, maxDepth)
                        });
                    }
                }
            }
            if (params.contains("ground-texture")) {
                model.textureClass = params.at("ground-texture").get<std::string>();
            }
        }

        for (auto& child : json_content.at("children")) {
            int childrenCount = 1; // Default value
            if (child.contains("params") && child.at("params").contains("quantity"))
                childrenCount = child.at("params").at("quantity").get<int>();

            for (int i = 0; i < childrenCount; i++) {
                if (child.contains("model-name") || child.contains("class")) {
                    model.modelChildren.push_back(BiomeModel::fromJson(child));
                } /*else if (child.contains("class")) {
                    model.instanceChildren.push_back(BiomeInstance::fromClass(child.at("class").get<std::string>()));
                }*/
            }
        }
    }
    return model;
}

std::shared_ptr<BiomeInstance> recursivelyCreateBiomeInstanceFromModel(BiomeModel model, Vector3 biomePosition, ShapeCurve area) {
    std::string biomeClass = model.modelName;
    // Should be able to retrieve the parameters of the biome...
    std::shared_ptr<BiomeInstance> instance = std::make_shared<BiomeInstance>(BiomeInstance::fromClass(biomeClass));
    instance->position = biomePosition;
    instance->area = area;
    instance->depthShape = model.depthShape;
    instance->textureClass = model.textureClass;
    auto children = model.modelChildren;

    Voronoi diagram(children.size(), area);
    if (instance->classname == "Mayotte") {
        std::ifstream file("C:/codes/Qt/generation-terrain-procedural/saved_maps/neighboring_constraints.json");
        std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
        nlohmann::json neighboring = nlohmann::json::parse(content);
        ConstraintsSolver solver;

        // Put the island in the middle
        diagram.pointset[0] = BSpline(diagram.pointset).center();

        float repulsionDistance = 100.f;
        float attractionDistance = 1.f;
        float neutralDistance = 50.f;
        for (size_t i = 0; i < diagram.pointset.size(); i++) {
            solver.addItem(new Vector3(diagram.pointset[i]));
        }
        std::vector<float> constraintsWeights(diagram.pointset.size());
        for (size_t i = 0; i < diagram.pointset.size(); i++) {
            if (neighboring.contains(model.modelChildren[i].modelName)) {
                auto biomeConstraints = neighboring.at(model.modelChildren[i].modelName);
                std::vector<std::string> cannot;
                std::vector<std::string> must;
                if (biomeConstraints.contains("cannot"))
                    cannot = biomeConstraints.at("cannot").get<std::vector<std::string>>();
                if (biomeConstraints.contains("must"))
                    must = biomeConstraints.at("must").get<std::vector<std::string>>();
                for (size_t j = 0; j < diagram.pointset.size(); j++) {
                    if (std::find(cannot.begin(), cannot.end(), model.modelChildren[j].modelName) != cannot.end()) {
                        constraintsWeights[i] += repulsionDistance;
                    }
                    else if (std::find(must.begin(), must.end(), model.modelChildren[j].modelName) != must.end()) {
                        constraintsWeights[i] += attractionDistance;
                    }
                    else {
                        constraintsWeights[i] += neutralDistance;
                    }
                }
            }
        }
        for (size_t i = 0; i < diagram.pointset.size(); i++) {
            if (neighboring.contains(model.modelChildren[i].modelName)) {
                auto biomeConstraints = neighboring.at(model.modelChildren[i].modelName);
                std::vector<std::string> cannot;
                std::vector<std::string> must;
                if (biomeConstraints.contains("cannot"))
                    cannot = biomeConstraints.at("cannot").get<std::vector<std::string>>();
                if (biomeConstraints.contains("must"))
                    must = biomeConstraints.at("must").get<std::vector<std::string>>();
                for (size_t j = 0; j < diagram.pointset.size(); j++) {
                    if (std::find(cannot.begin(), cannot.end(), model.modelChildren[j].modelName) != cannot.end()) {
                        solver.addDistanceConstraint(i, j, repulsionDistance / constraintsWeights[i]);
                    }
                    else if (std::find(must.begin(), must.end(), model.modelChildren[j].modelName) != must.end()) {
                        solver.addDistanceConstraint(i, j, attractionDistance / constraintsWeights[i]);
                    }
                    else {
                        solver.addDistanceConstraint(i, j, neutralDistance / constraintsWeights[i]);
                    }
                }
            }
        }
        auto newPositions = solver.solve(false, .1f, .1f);
        for (size_t i = 0; i < diagram.pointset.size(); i++) {
            diagram.pointset[i] = newPositions[i];
        }
        BSpline group(diagram.pointset);
        Vector3 AABBoxMin, AABBoxMax;
        std::tie(AABBoxMin, AABBoxMax) = group.AABBox();
        Vector3 containingBoxSize = group.containingBoxSize();
        for (size_t i = 0; i < diagram.pointset.size(); i++) {
            diagram.pointset[i] = ((diagram.pointset[i] - AABBoxMin) / containingBoxSize) * area.containingBoxSize();
        }
    }

    std::vector<BSpline> subarea_borders = diagram.solve();

    for (size_t i = 0; i < children.size() && i < diagram.pointset.size(); i++) {
//        std::cout << "Shape for " << children[i].modelName << " (son of " << model.modelName << ") :\n";
//        for (auto& p : subarea_borders[i].points)
//            std::cout << "- " << p << "\n";
//        std::cout << std::endl;
        std::shared_ptr<BiomeInstance> childBiome = recursivelyCreateBiomeInstanceFromModel(children[i], diagram.pointset[i], subarea_borders[i]);
        childBiome->parent = instance;
        instance->instances.push_back(childBiome);
    }
    return instance;
}
BiomeInstance BiomeModel::createInstance(Vector3 initialPosition, ShapeCurve initialArea)
{
    BiomeInstance instance = *recursivelyCreateBiomeInstanceFromModel(*this, initialPosition, initialArea);
    instance.completeIfNeeded();
    return instance;
}
