#include "EnvironmentalScene.h"



EnvironmentalScene::EnvironmentalScene()
{
    initFlow();
    this->terrainNormals = GridV3();
    this->materials = std::map<std::string, EnvMaterial>();
    this->availableObjects = std::map<std::string, EnvObject*>();
    this->instantiatedObjects = std::vector<EnvObject*>();
    this->flowImpactFactor = .9f;
    this->currentMaxID = -1;
    this->transformationRules = std::vector<MaterialsTransformation>();
    this->currentTime = 0;
    this->scenario = Scenario(this);

    this->allVectorProperties = std::map<std::string, GridV3>();
    this->allScalarProperties = std::map<std::string, GridF>();
}

GridV3 EnvironmentalScene::initFlow(bool force) {
    if (force || this->flowfield.empty()) {
        this->initialFlowfield = GridV3(100, 100, 1, Vector3(0, 0, 0));
        this->initialFlowfield.raiseErrorOnBadCoord = false;
        this->initialFlowfield.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;
        this->flowfield = this->initialFlowfield;
    }
    return this->flowfield;
}

void EnvironmentalScene::readEnvObjectsFile(std::string filename)
{
    std::ifstream file(filename);
    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    this->readEnvObjectsFileContent(content);
}

void EnvironmentalScene::readEnvObjectsFileContent(std::string content)
{
    auto json = nlohmann::json::parse(toLower(content));
    for (auto& obj : json) {
        std::string objName = obj["name"];
        if (startsWith(objName, "--")) continue; // Ignore some objects if the name starts with "--"
        if (!obj.contains("type")) {
            throw std::domain_error("No type given for Environmental Object defined as " + nlohmann::to_string(obj));
        }

        if (obj["type"] == "point")
            this->availableObjects[objName] = EnvPoint::fromJSON(obj);
        else if (obj["type"] == "curve")
            this->availableObjects[objName] = EnvCurve::fromJSON(obj);
        else if (obj["type"] == "area")
            this->availableObjects[objName] = EnvArea::fromJSON(obj);
        else {
            throw std::domain_error("Unrecognized type for Environmental Object defined as " + nlohmann::to_string(obj));
        }
    }


    for (auto& [name, obj] : this->availableObjects) {
        obj->fittingFunction = EnvObject::parseFittingFunction(obj->s_FittingFunction, obj->name, this);
        obj->fitnessFunction = EnvObject::parseFittingFunction(obj->s_FitnessFunction, obj->name, this);

        obj->snake.imageField = obj->fittingFunction;
        obj->snake.gradientField = gradientFromFieldFunction(obj->snake.imageField);
    }

    for (auto& obj : this->instantiatedObjects) {
        auto name = obj->name;
        obj->flowEffect = this->availableObjects[name]->flowEffect;
        //        obj->sandEffect = this->availableObjects[name]->sandEffect;
        obj->materialAbsorptionRate = this->availableObjects[name]->materialAbsorptionRate;
        obj->materialDepositionRate = this->availableObjects[name]->materialDepositionRate;
    }
    //    precomputeTerrainProperties(Heightmap());
}

void EnvironmentalScene::readEnvMaterialsFile(std::string filename)
{
    std::ifstream file(filename);
    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    this->readEnvMaterialsFileContent(content);
}

void EnvironmentalScene::readEnvMaterialsFileContent(std::string content)
{
    auto json = nlohmann::json::parse(toLower(content));
    for (auto& mat : json) {
        std::string matName = mat["name"];
        if (startsWith(matName, "--")) continue; // Ignore some materials if the name starts with "--"

        float diffusionSpeed = mat["diffusionspeed"];
        float waterTransport = mat["watertransport"];
        float mass = mat["mass"];
        float decay = mat["decay"];
        float virtualHeight = mat["virtualheight"];

        EnvMaterial material;
        material.name = matName;
        material.diffusionSpeed = diffusionSpeed;
        material.waterTransport = waterTransport;
        material.mass = mass;
        material.decay = decay;
        material.virtualHeight = virtualHeight;

        if (this->materials.count(matName) != 0) {
            material.currentState = this->materials[matName].currentState;
        } else {
            material.currentState = GridF(this->flowfield.getDimensions(), 0.f);
        }

        this->materials[matName] = material;
    }
}

void EnvironmentalScene::readEnvMaterialsTransformationsFile(std::string filename)
{
    std::ifstream file(filename);
    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    this->readEnvMaterialsTransformationsFileContent(content);
}

void EnvironmentalScene::readEnvMaterialsTransformationsFileContent(std::string content)
{
    std::vector<MaterialsTransformation> rules;
    //    std::string sline;
    auto lines = split(content, "\n");
    //    while (std::getline(content, sline)) {
    for (std::string sline : lines) {
        if (sline.empty() || sline[0] == '#') continue; // Comments with "#"
        std::istringstream line(sline);
        std::map<std::string, float> inputs, outputs;
        std::string value, word, operation;
        bool transformationValid = true;
        // Get inputs
        while (operation != "=") {
            line >> value;
            line >> word;
            line >> operation;
            inputs[word] = std::stof(value);
            transformationValid &= (this->materials.count(word) != 0);
        }
        // Get outputs
        while (true) {
            line >> value;
            line >> word;
            outputs[word] = std::stof(value);
            transformationValid &= (this->materials.count(word) != 0);
            if (!(line >> operation)) break;
        }
        if (transformationValid)
            rules.push_back({inputs, outputs});
        else {
            std::cerr << "Transformation not valid : " << sline << std::endl;
        }
    }
    this->transformationRules = rules;
}

void EnvironmentalScene::readScenarioFile(std::string filename)
{
    std::ifstream file(filename);
    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    this->readScenarioFileContent(content);
}

void EnvironmentalScene::readScenarioFileContent(std::string content)
{
    auto json = nlohmann::json::parse(toLower(content));

    scenario = Scenario(this);

    auto objects = json["objects"].get<std::map<std::string, nlohmann::json>>();
    for (auto [name, obj] : objects) {
        float proba = obj["proba"];
        int amount = (obj.contains("amount") ? obj["amount"].get<int>() : -1);

        scenario.addObject(name, proba, amount);
    }

    auto events = json["events"];
    for (auto& event : events) {
        std::string type = toLower(event["type"]);
        float startTime = event["start"];
        float endTime = event["end"];

        if (type == "waterlevel") {
            float amount = event["amount"];
            scenario.waterLevelEvents.push_back(WaterLevelEvent(amount, startTime, endTime));
        } else if (type == "storm") {
            Vector3 position = json_to_vec3<float>(event["position"]);
            Vector3 direction = json_to_vec3<float>(event["direction"]);
            float sigma = event["sigma"];
            scenario.stormEvents.push_back(StormEvent(position, direction, sigma, startTime, endTime));
        } else if (type == "subsidence") {
            Vector3 position = json_to_vec3<float>(event["position"]);
            float amount = event["amount"];
            float sigma = event["sigma"];
            scenario.subsidenceEvents.push_back(SubsidenceEvent(position, amount, sigma, startTime, endTime));
        } else if (type == "tectonic") {
            Vector3 direction = json_to_vec3<float>(event["direction"]);
            float sigma = event["sigma"];
            scenario.tectonicEvents.push_back(TectonicEvent(direction, sigma, startTime, endTime));
        } else {
            std::cerr << "The event " << type << " is not recognized..." << std::endl;
        }
    }

    auto parameters = json["simulation"];
    float duration = parameters["end"];
    float dt = parameters["dt"];
    float waterLevel = parameters["waterlevel"];

    scenario.duration = duration;
    scenario.dt = dt;
    scenario.waterLevel = waterLevel;

}


EnvObject *EnvironmentalScene::findClosest(std::string objectName, const Vector3 &pos)
{
    float minDist = std::numeric_limits<float>::max();
    EnvObject* bestElem = nullptr;
    for (auto& instance : this->instantiatedObjects) {
        if (instance->name != objectName) continue;
        float distance = instance->getSqrDistance(pos);
        if (distance < minDist) {
            minDist = distance;
            bestElem = instance;
        }
    }
    return bestElem;
}


EnvObject *EnvironmentalScene::instantiate(std::string objectName)
{
    if (this->availableObjects.count(objectName) == 0) {
        return nullptr;
    }
    this->currentMaxID++;
    auto object = this->availableObjects[objectName]->clone();
    object->ID = this->currentMaxID;
    object->scene = this;
    this->instantiatedObjects.push_back(object);
    return object;
}

void EnvironmentalScene::removeObject(EnvObject *obj)
{
    if (obj) {
        auto& list = this->instantiatedObjects;
        list.erase(std::find(list.begin(), list.end(), obj));
    }
}

void EnvironmentalScene::removeAllObjects()
{
    for (auto& object : this->instantiatedObjects) {
        delete object;
    }
    this->instantiatedObjects.clear();
}

bool EnvironmentalScene::applyEffects(const GridF& heights, const GridV3& userFlow)
{
    this->updateFlowfield(userFlow);
    return this->updateSedimentation(heights);
}

bool EnvironmentalScene::updateSedimentation(const GridF& heights)
{
    std::cout << "SEDIMENTATION" << std::endl;
    bool bigChangesInAtLeastOneMaterialDistribution = false;
    GridV3 heightsGradients = heights.gradient();
    auto smoothFluids = this->flowfield.meanSmooth(3, 3, 1, true);

    std::vector<std::string> names;
    /*for (auto& [name, material] : this->materials) {
        // TODO : SELECT ONLY AFFECTED MATERIALS
        // if (!material.isStable)
        names.push_back(name);
    }*/

    for (size_t i = 0; i < this->instantiatedObjects.size(); i++) {
        auto& object = this->instantiatedObjects[i];
        for (auto& [materialName, rate] : object->materialAbsorptionRate) {
            if (rate > 0 && std::find(names.begin(), names.end(), materialName) != names.end()) {
                names.push_back(materialName);
            }
        }
        for (auto& [materialName, rate] : object->materialDepositionRate) {
            if (rate > 0 && std::find(names.begin(), names.end(), materialName) != names.end()) {
                names.push_back(materialName);
            }
        }
    }
#pragma omp parallel for
    for(int i = 0; i < names.size(); i++) {
        auto& material = materials[names[i]];
        // if (material.isStable) {
        // std::cout << material.name << " is stable" << std::endl;
        // continue;
        // } else {
        // std::cout << material.name << " NOT stable" << std::endl;
        // }

        // bool needToBeUpdated = false;
        // float startingAmount = material.currentState.sum();
        auto startState = material.currentState;
        for (size_t i = 0; i < this->instantiatedObjects.size(); i++) {
            auto& object = this->instantiatedObjects[i];
            if (object->materialAbsorptionRate.count(material.name) != 0 && object->materialAbsorptionRate[material.name] != 0) {
                // #pragma omp critical
                object->applyAbsorption(material);
                // needToBeUpdated = true;
            }
        }

        for (size_t i = 0; i < this->instantiatedObjects.size(); i++) {
            auto& object = this->instantiatedObjects[i];
            if (object->materialDepositionRate.count(material.name) != 0 && object->materialDepositionRate[material.name] != 0) {
                // #pragma omp critical
                object->applyDeposition(material);
                // needToBeUpdated = true;
            }
        }

        // if (needToBeUpdated) {
        material.update(smoothFluids, heightsGradients, this->scenario.dt);
        // material.currentState *= material.decay;

        // float endingAmount = material.currentState.sum();
        // float diff = (material.currentState - startState).abs().sum();

        // #pragma omp critical
        // {
        //     std::cout << material.name << ": " << diff << std::endl;
        // }

        // if (diff > 1e-3) {
        // #pragma omp critical
        // {
        // bigChangesInAtLeastOneMaterialDistribution = true;
        // }
        // } else {
        // std::cout << material.name << " diff : " << std::abs(endingAmount - startingAmount) << std::endl;
        // material.isStable = true;
        // }
        // } else {
        material.isStable = true;
        // }
    }
    return bigChangesInAtLeastOneMaterialDistribution;
}

std::vector<std::string> EnvironmentalScene::updateSedimentationKnowingFluidsAndGradients([[maybe_unused]] const GridF& heights, const GridV3& heightsGradients, const GridV3& smoothFluids, std::vector<std::string> unstableMaterials)
{
    std::cout << "SED_FLUID" << std::endl;
    std::vector<std::string> stillUnstable;
    std::vector<std::string> names = unstableMaterials;

#pragma omp parallel for
    for(int i = 0; i < names.size(); i++) {
        auto& material = materials[names[i]];

        bool needToBeUpdated = false;
        auto startState = material.currentState;
        for (size_t i = 0; i < this->instantiatedObjects.size(); i++) {
            auto& object = this->instantiatedObjects[i];
            if (object->materialAbsorptionRate.count(material.name) != 0 && object->materialAbsorptionRate[material.name] != 0) {
                object->applyAbsorption(material);
                needToBeUpdated = true;
            }
        }

        for (size_t i = 0; i < this->instantiatedObjects.size(); i++) {
            auto& object = this->instantiatedObjects[i];
            if (object->materialDepositionRate.count(material.name) != 0 && object->materialDepositionRate[material.name] != 0) {
                object->applyDeposition(material);
                needToBeUpdated = true;
            }
        }

        // if (needToBeUpdated) {
        material.update(smoothFluids, heightsGradients, this->scenario.dt);
        float diff = (material.currentState - startState).abs().sum();
        if (diff > 1e-3) {
#pragma omp critical
            {
                stillUnstable.push_back(material.name);
            }
        }
        // }
    }
    return stillUnstable;
}

void EnvironmentalScene::stabilizeMaterials(const GridF &heights, int maxIterations)
{
    GridV3 heightsGradients = heights.gradient();
    auto smoothFluids = this->flowfield.meanSmooth(3, 3, 1, true);
    std::vector<std::string> unstableMaterials;
    for (auto& [name, material] : this->materials) {
        unstableMaterials.push_back(name);
    }

    for (int iteration = 0; iteration < maxIterations; iteration++) {
        unstableMaterials = this->updateSedimentationKnowingFluidsAndGradients(heights, heightsGradients, smoothFluids, unstableMaterials);
        if (unstableMaterials.empty()) {
            break;
        }
    }
}

void EnvironmentalScene::applyMaterialsTransformations()
{
    displayProcessTime("Filling compact materials... ", [&]() {
            std::set<std::string> neededMaterials;
            for (size_t iRule = 0; iRule < transformationRules.size(); iRule++) {
                auto [input, output] = transformationRules[iRule];
                for (auto [inMaterial, inDose] : input) {
                    neededMaterials.insert(inMaterial);
                }
                for (auto [outMaterial, outDose] : output) {
                    neededMaterials.insert(outMaterial);
                }
            }
            std::map<std::string, float> initialState; // Loop the map creation only once
            for (const auto& matName : neededMaterials)
                initialState.insert({matName, 0.f});
            Matrix3<std::map<std::string, float>> allMaterials(this->flowfield.getDimensions(), initialState);
            allMaterials.iterateParallel([&] (size_t i) {
                for (const auto& [matName, amount] : allMaterials[i]) {
                    allMaterials[i][matName] = this->materials[matName].currentState[i];
                }

                for (size_t iRule = 0; iRule < transformationRules.size(); iRule++) {
                    const auto& [input, output] = transformationRules[iRule];
                    float maxTransform = 10000.f;
                    for (const auto& [inMaterial, inDose] : input) {
                        float inAmount = allMaterials[i][inMaterial];
                        float transformVal = inAmount / inDose;
                        maxTransform = std::min(maxTransform, transformVal);
                    }
                    if (maxTransform > 1e-3) {
                        for (const auto& [inMaterial, inDose] : input) {
                            allMaterials[i][inMaterial] -= inDose * maxTransform;
                        }
                        for (const auto& [outMaterial, outDose] : output) {
                            allMaterials[i][outMaterial] += outDose * maxTransform;
                        }
                    }
                }
            });

            for (auto& matName : neededMaterials) {
                auto& mat = this->materials[matName];
                mat.currentState.iterateParallel([&](size_t i) {
                    mat.currentState[i] = allMaterials[i][matName];
                });
            }
        }, false);
}

void EnvironmentalScene::updateFlowfield(const GridV3 &userFlow)
{
    this->flowfield = this->initialFlowfield;
    if (!userFlow.empty())
        this->flowfield += userFlow;
    for (int i = 0; i < this->instantiatedObjects.size(); i++) {
        auto& object = this->instantiatedObjects[i];
        auto [flow, occupancy] = object->computeFlowModification();
        this->flowfield = flow;
    }
    this->flowfield = this->flowfield.meanSmooth(3, 3, 1, true);
}

void EnvironmentalScene::beImpactedByEvents()
{
    for (auto& obj : this->instantiatedObjects) {
        obj->age += 1.f;
    }
}


void EnvironmentalScene::precomputeTerrainProperties(const GridF& heightmap, float waterLevel, float maxHeight)
{

    displayProcessTime("Computing terrain properties... ", [&]() {
        Vector3 terrainDimensions = heightmap.getDimensions();
        GridF initialScalarPropertyMap(terrainDimensions, 0.f);
        GridV3 initialVectorPropertyMap(terrainDimensions, Vector3::invalid());

        // Initialize the maps
        for (auto& [name, obj] : this->availableObjects) {
            this->allVectorProperties[name] = initialVectorPropertyMap;
            this->allVectorProperties[name + ".center"] = initialVectorPropertyMap;
            this->allVectorProperties[name + ".start"] = initialVectorPropertyMap;
            this->allVectorProperties[name + ".end"] = initialVectorPropertyMap;
            this->allVectorProperties[name + ".normal"] = initialVectorPropertyMap;
            this->allVectorProperties[name + ".dir"] = initialVectorPropertyMap;
            this->allScalarProperties[name + ".inside"] = initialScalarPropertyMap;
            this->allScalarProperties[name + ".curvature"] = initialScalarPropertyMap;
        }
        this->allVectorProperties["current"] = initialVectorPropertyMap;
        this->allVectorProperties["current.dir"] = initialVectorPropertyMap;
        this->allScalarProperties["current.vel"] = initialScalarPropertyMap;
        this->allVectorProperties["current.gradient"] = initialVectorPropertyMap;

        this->allScalarProperties["depth"] = initialScalarPropertyMap;
        this->allVectorProperties["depth.gradient"] = initialVectorPropertyMap;

        this->allScalarProperties["fracture"] = initialScalarPropertyMap;
        this->allVectorProperties["fracture.gradient"] = initialVectorPropertyMap;

        for (const auto& [matName, material] : this->materials) {
            this->allScalarProperties[matName] = initialScalarPropertyMap;
            this->allVectorProperties[matName + ".gradient"] = initialVectorPropertyMap;
        }


        // Evaluate at each point
        for (auto& [name, obj] : this->availableObjects) {
            displayProcessTime("Computing properties for " + name + "... ", [&]() {
                    this->recomputeTerrainPropertiesForObject(name);
                }, false);
        }
        this->recomputeFlowAndSandProperties(heightmap, waterLevel, maxHeight);
    });
}

void EnvironmentalScene::recomputeTerrainPropertiesForObject(std::string objectName)
{
    auto name = objectName;
    this->flowfield.iterateParallel([&](const Vector3i& pos) {
        //        auto [distance, object] = this->getSqrDistanceTo(name, pos);
        EnvObject* object = this->findClosest(objectName, pos);
        if (object == nullptr) {
            this->allVectorProperties[name](pos) = Vector3::invalid();
            this->allVectorProperties[name + ".center"](pos) = Vector3::invalid();
            this->allVectorProperties[name + ".start"](pos) = Vector3::invalid();
            this->allVectorProperties[name + ".end"](pos) = Vector3::invalid();
            this->allScalarProperties[name + ".inside"](pos) = 0.f;
            this->allVectorProperties[name + ".normal"](pos) = Vector3::invalid();
            this->allVectorProperties[name + ".dir"](pos) = Vector3::invalid();
            this->allScalarProperties[name + ".curvature"](pos) = 0.f;
        } else {
            auto allProperties = object->getAllProperties(pos);
            this->allVectorProperties[name](pos) = allProperties["default"];
            this->allVectorProperties[name + ".center"](pos) = allProperties["center"];
            this->allVectorProperties[name + ".start"](pos) = allProperties["start"];
            this->allVectorProperties[name + ".end"](pos) = allProperties["end"];
            this->allScalarProperties[name + ".inside"](pos) = (allProperties["inside"].isValid() ? 1.f : 0.f);
            this->allVectorProperties[name + ".normal"](pos) = allProperties["normal"];
            this->allVectorProperties[name + ".dir"](pos) = allProperties["dir"];
            this->allScalarProperties[name + ".curvature"](pos) = (allProperties["curvature"].x() < 1e5 ? allProperties["curvature"].x() : -1.f);
        }
    });
}

void EnvironmentalScene::recomputeFlowAndSandProperties(const GridF& heightmap, float waterLevel, float maxHeight)
{
    this->recomputeFlow();
    for (auto& [matName, material] : this->materials) {
        this->allScalarProperties[matName] = material.currentState;
        this->allVectorProperties[matName + ".gradient"] = material.currentState.gradient();
    }
    this->allScalarProperties["depth"] = ((waterLevel * maxHeight) - heightmap.meanSmooth(3, 3, 1));
    // this->allScalarProperties["depth"] = ((waterLevel * maxHeight) - heightmap.gaussianSmooth(1.f, true));
    this->allVectorProperties["depth.gradient"] = this->allScalarProperties["depth"].gradient();

    this->allScalarProperties["fracture"] = this->scenario.computeTectonic(this->allScalarProperties["fracture"].getDimensions());
    this->allVectorProperties["fracture.gradient"] = this->allScalarProperties["fracture"].gradient();
}

void EnvironmentalScene::recomputeFlow()
{
    this->flowfield.iterateParallel([&](const Vector3i& pos) {
        Vector3 waterFlow = this->flowfield(pos);
        this->allVectorProperties["current"](pos) = waterFlow;
        this->allVectorProperties["current.dir"](pos) = waterFlow.normalized();
        this->allScalarProperties["current.vel"](pos) = waterFlow.length();
    });
    this->allVectorProperties["current.gradient"] = this->allScalarProperties["current.vel"].gradient();
}

GridF EnvironmentalScene::getHeightmap(const GridF& initialHeightmap, float absoluteWaterLevel, float flowErosionFactor, bool displayGrooves)
{
    GridF subsidedHeightmap = GridF(initialHeightmap.getDimensions());
    GridF subsidenceFactor = this->scenario.computeSubsidence(initialHeightmap.getDimensions());
    subsidedHeightmap = initialHeightmap * subsidenceFactor;

    GridF groundConstraintedHeights = GridF(subsidedHeightmap.getDimensions()); // Heightmaps from the ground
    GridF waterConstraintedHeights = GridF(subsidedHeightmap.getDimensions(), -100000.f); // Heightmaps from the water level
    GridF surfaceHeights = GridF(subsidedHeightmap.getDimensions());
    for (auto& obj : this->instantiatedObjects) {
        if (auto patch = dynamic_cast<ImplicitPrimitive*>(obj->_patch)) {
            GridF grid = GridF(subsidedHeightmap.getDimensions(), 0.f);
            grid = grid.paste(obj->createHeightfield() * obj->computeGrowingState2(), patch->position.xy());
            if (flowErosionFactor != 0 && this->materials.count(toLower(stringFromMaterial(obj->material)))) {
                grid = grid.warpWith(this->flowfield * flowErosionFactor * this->materials[toLower(stringFromMaterial(obj->material))].waterTransport, 10);
            }
            if (obj->heightFrom == EnvObject::HeightmapFrom::SURFACE) {
                surfaceHeights = (surfaceHeights + grid * (isIn(obj->material, LayerBasedGrid::invisibleLayers) ? -1.f : 1.f)).max(-15.f);
            } else if (obj->heightFrom == EnvObject::HeightmapFrom::GROUND) {
                groundConstraintedHeights = groundConstraintedHeights.max(grid * subsidenceFactor, Vector3());
            } else if (obj->heightFrom == EnvObject::HeightmapFrom::WATER) {
                grid.iterateParallel([&] (size_t i) {
                    grid[i] = (std::abs(grid[i]) < 1e-4 ? -10000.f : grid[i]);
                });
                // std::cout << "Max height for " << obj->name << ": " << grid.max() << " while height = " << obj->height << "(grow = " <<  obj->computeGrowingState2() << ")" << std::endl;
                waterConstraintedHeights = waterConstraintedHeights.max((grid - (obj->height)) - (obj->name == "lagoon" || obj->name == "smalllagoon" ? 3.f : 1.f), Vector3()); // Not sure why I need to multiply by 2.0, but otherwise, maxHeight is heigher than obj->height...
            }
        }
    }
    // if (flowErosionFactor != 0) {
    // groundConstraintedHeights = groundConstraintedHeights.warpWith(this->flowfield * flowErosionFactor, 10);
    // waterConstraintedHeights = waterConstraintedHeights.warpWith(this->flowfield * flowErosionFactor, 10);
    // surfaceHeights = surfaceHeights.warpWith(this->flowfield * flowErosionFactor, 10);
    // }
    // Dirty, remove when you understand why lagoon get over the water...
    waterConstraintedHeights.iterateParallel([&](size_t i) {
        waterConstraintedHeights[i] = std::min(waterConstraintedHeights[i] + absoluteWaterLevel, absoluteWaterLevel - 1.f);
    });
    waterConstraintedHeights = waterConstraintedHeights.meanSmooth(3, 3, 1);
    // waterConstraintedHeights = waterConstraintedHeights.gaussianSmooth(1.f, true);

    bool modificationsAppliedToSurface = false;
    for (auto& obj : this->instantiatedObjects) {
        if (displayGrooves) {
            if (endsWith(toLower(obj->name), "reef")) {
                auto objAsEnvCurve = dynamic_cast<EnvCurve*>(obj);
                BSpline path = objAsEnvCurve->curve;
                float nbGrooves = path.length() / 10.f;
                float sigma = objAsEnvCurve->width;
                surfaceHeights.iterateParallel([&](const Vector3i& pos) {
                    float closestT = path.estimateClosestTime(pos);
                    float closestGrooveStartT = float(int(closestT * nbGrooves)) / nbGrooves;
                    auto [closestPoint, direction, normal] = path.pointAndDerivativeAndSecondDerivative(closestT);
                    auto closestGrooveStartPoint = path.getPoint(closestGrooveStartT);
                    if (direction.norm2() == 0) return;
                    direction.normalize();
                    auto fakeNormal = direction.rotated90XY(); // (normal.norm2() > 0 ? normal.normalize() : direction.rotated90XY());
                    Vector3 newSpace = Vector3(pos - closestGrooveStartPoint).changeBasis(direction, fakeNormal, Vector3(0, 0, 1)); //.rotated(Vector3(0, 0, random_gen::generate_perlin(closestT * 500.f) * 0.2f));
                    float sizeX = 1.f/(nbGrooves * .5f), sizeY = 1.f/(sigma * 1.f);
                    // float initialDistance = std::clamp(1.f - (pos - closestPoint).norm() / sigma, 0.f, 1.f);
                    float grooves = std::max(0.f, 1.f - (sizeX * std::abs(newSpace.x() - 1.f/sizeX) + std::pow(sizeY * newSpace.y(), 2.f)));
                    // return std::max(grooves, initialDistance);
                    const Vector3& flow = this->flowfield(pos);
                    surfaceHeights(pos) += 2.f * grooves * std::max(abs(flow.dot(fakeNormal)), 0.f);
                });
                modificationsAppliedToSurface = true;
            }
        }
    }

    if (modificationsAppliedToSurface) {
        surfaceHeights = surfaceHeights.meanSmooth(3, 3, 1);
        // surfaceHeights = surfaceHeights.gaussianSmooth(1.f, true, true);
    }

    subsidedHeightmap = GridF::max(GridF::max(subsidedHeightmap, groundConstraintedHeights), waterConstraintedHeights).meanSmooth(5, 5, 1);
    subsidedHeightmap = (subsidedHeightmap.max(-15.f) + surfaceHeights).meanSmooth(3, 3, 1).max(-15.f);
    // subsidedHeightmap = GridF::max(GridF::max(subsidedHeightmap, groundConstraintedHeights), waterConstraintedHeights).gaussianSmooth(2.f, true, true);
    // subsidedHeightmap = (subsidedHeightmap.max(-15.f) + surfaceHeights).gaussianSmooth(1.f, true, true).max(-15.f);
    return subsidedHeightmap;
}

void EnvironmentalScene::reset()
{
    std::set<std::string> destroyedObjects;
    for (auto& obj : this->instantiatedObjects) {
        destroyedObjects.insert(obj->name);
        delete obj;
    }
    this->instantiatedObjects.clear();
    for (auto name : destroyedObjects) {
        this->recomputeTerrainPropertiesForObject(name);
    }
    initFlow(true);
    for (auto& [matName, mat] : materials) {
        mat.currentState.reset();
    }

}

/*
#include "Utils/Delaunay.h"
GraphObj EnvironmentalScene::sceneToGraph()
{
    GraphObj graph;
    std::vector<GraphNodeObj*> nodes(this->instantiatedObjects.size());
    std::vector<Vector3> positions(nodes.size());

    for (int i = 0; i < nodes.size(); i++) {
        auto& obj = this->instantiatedObjects[i];
        positions[i] = dynamic_cast<EnvPoint*>(obj)->position;
        //        nodes[i] = graph.addNode(new GraphNodeObj(obj, positions[i], i));
    }

    graph = Delaunay().fromVoronoi(Voronoi(positions)).graph.cast<EnvObject*>();

    return graph;
}
*/
