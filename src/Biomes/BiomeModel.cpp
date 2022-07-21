#include "BiomeModel.h"
#include "Utils/AdjencySolver.h"
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
    instance->area = area.removeDuplicates();
    instance->depthShape = model.depthShape;
    instance->textureClass = model.textureClass;
    auto children = model.modelChildren;

    Voronoi diagram(children.size(), area);
    std::vector<BSpline> subarea_borders = diagram.solve(100); // Add some relaxations to be a little bit more uniform

    std::vector<std::string> allChildrenClassnames;
    for (size_t i = 0; i < children.size() && i < diagram.pointset.size(); i++) {
        allChildrenClassnames.push_back(children[i].modelName);
    }
#ifdef linux
    std::ifstream file("/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/saved_maps/neighboring_constraints.json");
#else
    std::ifstream file("C:/codes/Qt/generation-terrain-procedural/saved_maps/neighboring_constraints.json");
#endif
    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    nlohmann::json neighboring = nlohmann::json::parse(content);
    std::vector<float> constraintsWeights(diagram.pointset.size());
    std::vector<std::pair<std::string, std::string>> obligationNeighboring;
    std::vector<std::pair<std::string, std::string>> forbiddenNeighboring;
    for (nlohmann::json::iterator it = neighboring.begin(); it != neighboring.end(); ++it) {
        std::string biome1 = it.key();
        std::vector<std::string> obligations = it.value().at("must").get<std::vector<std::string>>();
        std::vector<std::string> forbidding = it.value().at("cannot").get<std::vector<std::string>>();
        for (auto& biome : obligations) {
            obligationNeighboring.push_back(std::make_pair(biome1, biome));
            obligationNeighboring.push_back(std::make_pair(biome, biome1));
        }
        for (auto& biome : forbidding) {
            forbiddenNeighboring.push_back(std::make_pair(biome1, biome));
            forbiddenNeighboring.push_back(std::make_pair(biome, biome1));
        }
    }
    std::vector<int> newChildrenOrder = std::vector<int>(children.size());
    for (size_t i = 0; i < newChildrenOrder.size(); i++) newChildrenOrder[i] = i;
    AdjencySolver solver;
//    auto start = std::chrono::system_clock::now();
    newChildrenOrder = solver.solve(diagram, allChildrenClassnames /*vectorMerge(allChildrenClassnames, allChildrenClassnames)*/, forbiddenNeighboring);
    // Because we may have used multiple of the data, we need to clamp the new indices
    for (auto& ind : newChildrenOrder) ind = ind % allChildrenClassnames.size();
//    auto end = std::chrono::system_clock::now();
//    std::cout << "Result in " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "s" << std::endl;

    for (size_t i = 0; i < children.size() && i < diagram.pointset.size(); i++) {
        int index = newChildrenOrder[i];
        std::shared_ptr<BiomeInstance> childBiome = recursivelyCreateBiomeInstanceFromModel(children[index], diagram.pointset[i], subarea_borders[i]);
        childBiome->parent = instance;
        instance->instances.push_back(childBiome);
    }
    return instance;
}
std::shared_ptr<BiomeInstance> BiomeModel::createInstance(Vector3 initialPosition, ShapeCurve initialArea)
{
    std::shared_ptr<BiomeInstance> instance = recursivelyCreateBiomeInstanceFromModel(*this, initialPosition, initialArea);
    instance->completeIfNeeded();
    return instance;
}
