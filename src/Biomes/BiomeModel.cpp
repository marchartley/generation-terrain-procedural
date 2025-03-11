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
                model.minDepth = depthInfo.at("min").get<float>();
                model.maxDepth = depthInfo.at("max").get<float>();
                if (depthInfo.contains("slope")) {
                    std::string slopeType = depthInfo.at("slope").get<std::string>();
                    if (slopeType == "mountain") {
                        model.slopeType = MOUNTAIN;
                    } else if (slopeType == "curve") {
                        model.slopeType = CURVE;
                    } else if (slopeType == "flat") {
                        model.slopeType = FLAT;
                    } else if (slopeType == "linear") {
                        model.slopeType = LINEAR_SLOPE;
                    } else  {
                        model.slopeType = UNDEFINED;
                    }
                    /*
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
                    }*/
                }
            }
            if (params.contains("ground-texture")) {
                model.textureClass = params.at("ground-texture").get<std::string>();
            }
        }

        for (auto& child : json_content.at("children")) {
            Probability probaAppearence = {1.0, 0.0}; // default value
            Probability probaSize;
            Probability probaQuantity = {1.0, 0.0}; // default value
            int priority = 0;
//            float childrenCount = 1; // Default value
            if (child.contains("params")) {
                if (child.at("params").contains("quantity")) {
                    probaQuantity = Probability::fromJson(child.at("params").at("quantity"));
                }
                if (child.at("params").contains("size")) {
                    probaSize = Probability::fromJson(child.at("params").at("size"));
                }
                if (child.at("params").contains("probability-appearence")) {
                    probaAppearence = Probability::fromJson(child.at("params").at("probability-appearence"));
                }
                if (child.at("params").contains("priority-offset")) {
                    priority = child.at("params").at("priority-offset").get<int>();
                }
            }
            std::shared_ptr<BiomeModel> childModel = std::make_shared<BiomeModel>(BiomeModel::fromJson(child));
            childModel->idealSize = probaSize;
            model.modelChildren.push_back({childModel, probaQuantity, probaAppearence, probaSize, priority});
//            for (int i = 0; i < childrenCount; i++) {
//                if (child.contains("model-name") || child.contains("class")) {
//                    model.modelChildren.push_back(BiomeModel::fromJson(child));
//                } /*else if (child.contains("class")) {
//                    model.instanceChildren.push_back(BiomeInstance::fromClass(child.at("class").get<std::string>()));
//                }*/
//            }
        }
    }
    return model;
}

nlohmann::json BiomeModel::toJson()
{
    nlohmann::json parameters;
    parameters["ground-texture"] = this->textureClass;

    if (this->slopeType != UNDEFINED) {
        std::string slope = "";
        if (slopeType == MOUNTAIN)     slope = "mountain";
        if (slopeType == CURVE)        slope = "curve";
        if (slopeType == FLAT)         slope = "flat";
        if (slopeType == LINEAR_SLOPE) slope = "linear";
        parameters["depth"] = nlohmann::json({{"min", this->minDepth}, {"max", this->maxDepth}, {"slope", slope}});
    }

    nlohmann::json children;
    for (auto& biomeChild : this->modelChildren) {
        nlohmann::json childJson = biomeChild.model->toJson();
        childJson["params"]["quantity"] = biomeChild.probaQuantity.toJson();
        childJson["params"]["size"] = biomeChild.probaSize.toJson();
        childJson["params"]["probability-appearence"] = biomeChild.probaAppearence.toJson();
        childJson["params"]["priority-offset"] = biomeChild.priorityOffset;
        children.push_back(childJson);
    }

    nlohmann::json content;
    content["model-name"] = this->modelName;
    content["params"] = parameters;
    if (!children.empty())
        content["children"] = children;

    return content;
}

std::shared_ptr<BiomeInstance> recursivelyCreateBiomeInstanceFromModel(std::shared_ptr<BiomeModel> model, const Vector3& biomePosition, ShapeCurve area) {
    std::string biomeClass = model->modelName;
    // Should be able to retrieve the parameters of the biome...
    std::shared_ptr<BiomeInstance> instance = std::make_shared<BiomeInstance>(BiomeInstance::fromClass(biomeClass));
    instance->model = model;
    instance->position = biomePosition;
    instance->priorityOffset = model->priortyOffset;
    instance->area = area.removeDuplicates();
    instance->depthShape = BSpline({Vector3(0, model->minDepth), Vector3(1, model->maxDepth)});
    instance->textureClass = model->textureClass;
    instance->idealSize = (model->idealSize.mean > 0 ? model->idealSize.randomValue() : -1);
    std::vector<std::shared_ptr<BiomeModel>> children;

    int augmentation = 1;
    for (size_t i = 0; i < model->modelChildren.size() * (augmentation + 1); i++) {
        auto& child = model->modelChildren[i % model->modelChildren.size()];
        for (int i = 0; i < child.probaQuantity.randomValue() * child.probaAppearence.randomValue(); i++) {
            children.push_back(child.model);
        }
    }

    Voronoi diagram(children.size(), area);
    std::vector<ShapeCurve> subarea_borders = diagram.solve(true, 100); // Add some relaxations to be a little bit more uniform

    std::vector<std::string> allChildrenClassnames;
    for (size_t i = 0; i < children.size() && i < diagram.pointset.size(); i++) {
            allChildrenClassnames.push_back(children[i]->modelName);
    }
    std::ifstream file("saved_maps/neighboring_constraints.json");
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

    newChildrenOrder = solver.solveGeomToTopo(diagram, allChildrenClassnames, forbiddenNeighboring);

//    std::vector<BiomeInstance> instances

    // Because we may have used multiple of the data, we need to clamp the new indices
//    for (auto& ind : newChildrenOrder)
//        ind = ind % allChildrenClassnames.size();

    std::vector<std::string> finalClasses(newChildrenOrder.size());
    for (size_t i = 0; i < finalClasses.size(); i++) {
        finalClasses[i] = allChildrenClassnames[newChildrenOrder[i]];
    }
    std::vector<int> mergingQueue = {0};

    std::vector<Vector3> newPointset = diagram.pointset;
    std::vector<ShapeCurve> newSubareas(subarea_borders.size());

    std::vector<std::vector<int>> groups;
    std::set<int> visitedNodes;
    std::vector<int> unvisitedNodes = newChildrenOrder;
    std::vector<int> queue;
    if (newChildrenOrder.size() > 1) { // Need at least 2 elems to have neighbors
        while (!unvisitedNodes.empty()) {
            groups.push_back({});
            queue = {unvisitedNodes.front()};
            while (!queue.empty()) {
                int currentIndex = queue.front();
                queue.erase(queue.begin());
                visitedNodes.insert(currentIndex);
                unvisitedNodes.erase(std::find(unvisitedNodes.begin(), unvisitedNodes.end(), currentIndex));
                groups[groups.size() - 1].push_back(currentIndex);
                auto currentClass = finalClasses[currentIndex];

                newSubareas[groups[groups.size() - 1][0]] = newSubareas[groups[groups.size() - 1][0]].merge(subarea_borders[currentIndex]);

                auto neighbors = diagram.neighbors[currentIndex];

                for (size_t iNeighbor = 0; iNeighbor < neighbors.size(); iNeighbor++) {
                    auto neighborClass = finalClasses[neighbors[iNeighbor]];
                    if (!isIn(neighbors[iNeighbor], queue) && visitedNodes.find(neighbors[iNeighbor]) == visitedNodes.end() && currentClass == neighborClass) {
                        queue.push_back(neighbors[iNeighbor]);
                        newPointset[neighbors[iNeighbor]].setValid(false);
                    }
                }
            }
        }
    }
    std::vector<std::set<int>> merged(newChildrenOrder.size());
    for (size_t iGrp = 0; iGrp < groups.size(); iGrp++) {
        merged[groups[iGrp][0]] = convertVectorToSet(groups[iGrp]);
    }
/*
    mergingQueue = newChildrenOrder;
    std::vector<std::set<int>> merged(mergingQueue.size());
    for (size_t i = 0; i < merged.size(); i++)
        merged[i] = {int(i)};
    if (newChildrenOrder.size() > 1) { // Need at least 2 elems to have neighbors
        while (!mergingQueue.empty()) {
            int currentIndex = mergingQueue.front();
            mergingQueue.erase(mergingQueue.begin());
            auto currentClass = finalClasses[currentIndex];
            auto neighbors = diagram.neighbors[currentIndex];

            for (size_t iNeighbor = 0; iNeighbor < neighbors.size(); iNeighbor++) {
                auto neighborClass = finalClasses[neighbors[iNeighbor]];
                if (currentClass == neighborClass) {
                    merged[currentIndex].insert(merged[neighbors[iNeighbor]].begin(), merged[neighbors[iNeighbor]].end());
                }
            }
        }
    }*/
    /*std::vector<Vector3> newPointset = diagram.pointset;
    std::vector<ShapeCurve> newSubareas(subarea_borders.size());
    for (size_t i = 0; i < merged.size(); i++) {
        if (merged[i].size() == 0) {
            newPointset[i].setValid(false);
            continue;
        }
        newSubareas[i] = ShapeCurve(subarea_borders[i]).removeDuplicates();
        for (auto neighbor: merged[i]) {
            if (neighbor != int(i)) {
                merged[neighbor].clear();
                newSubareas[i] = newSubareas[i].merge(ShapeCurve(subarea_borders[neighbor]));
            }
        }

    }*/
    for (int i = newPointset.size() - 1; i >= 0; i--) {
        if (!newPointset[i].isValid()) {
            newPointset.erase(newPointset.begin() + i);
            newSubareas.erase(newSubareas.begin() + i);
            newChildrenOrder.erase(newChildrenOrder.begin() + i);
        }
    }


    for (size_t i = 0; i < newChildrenOrder.size(); i++) {
        int index = newChildrenOrder[i];
        std::shared_ptr<BiomeInstance> childBiome = recursivelyCreateBiomeInstanceFromModel(children[index], newPointset[i], newSubareas[i]); // diagram.pointset[i], subarea_borders[i]);
        childBiome->parent = instance;
        instance->instances.push_back(childBiome);
    }
    return instance;
}
std::shared_ptr<BiomeInstance> BiomeModel::createInstance(const Vector3& initialPosition, ShapeCurve initialArea)
{
    std::shared_ptr<BiomeInstance> instance = recursivelyCreateBiomeInstanceFromModel(std::shared_ptr<BiomeModel>(this), initialPosition, initialArea);
    instance->completeIfNeeded();
    return instance;
}

std::shared_ptr<BiomeModel> BiomeModel::clone()
{
    std::shared_ptr<BiomeModel> model = std::make_shared<BiomeModel>(*this); //->shared_from_this();
    std::vector<BiomeModelChild> deepCopiedChildren(model->modelChildren.size());
    for (size_t i = 0; i < model->modelChildren.size(); i++) {
        deepCopiedChildren[i] = model->modelChildren[i];
        deepCopiedChildren[i].model = model->modelChildren[i].model->clone();
    }
    model->modelChildren = deepCopiedChildren;
    return model;
}
