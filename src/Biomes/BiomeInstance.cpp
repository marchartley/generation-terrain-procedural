#include "BiomeInstance.h"
#include "Utils/Voronoi.h"

std::map<int, std::shared_ptr<BiomeInstance>> BiomeInstance::instancedBiomes;

BiomeInstance::BiomeInstance()
{
    int ID = 0;
    for (auto& inst : BiomeInstance::instancedBiomes)
        ID = std::max(ID, inst.first + 1);
}

BiomeInstance BiomeInstance::fromClass(std::string className)
{
    BiomeInstance biome;
    biome.classname = className;
    return biome;
}

int BiomeInstance::getLevel(bool ignorePriorityOffset)
{
    if (!this->parent) {
        return 0 + (ignorePriorityOffset ? 0 : priorityOffset);
    } else {
        return this->parent->getLevel(true) + 1 + (ignorePriorityOffset ? 0 : priorityOffset);
    }
}

void BiomeInstance::completeIfNeeded()
{
    if (classname == "arche")
        completeArch();
    else if (classname == "tranchee" || classname == "passe-corail")
        completeTrench();
    else if (classname == "area")
        completeArea();

    for (auto& child : instances)
        child->completeIfNeeded();
}

std::shared_ptr<BiomeInstance> BiomeInstance::clone(ShapeCurve newArea, Vector3 newPosition)
{
    std::shared_ptr<BiomeInstance> cloneBiome = std::make_shared<BiomeInstance>(*this);
    BiomeInstance::registerBiomeInstance(cloneBiome);
    cloneBiome->area = newArea;
    if (!newPosition.isValid())
        newPosition = newArea.center();
    cloneBiome->position = newPosition;
//    std::cout << "Previous area : " << this->area << " - new area : " << cloneBiome->area << std::endl;
//    std::cout << "Previous instance size : " << this->instances.size() << " - new size : ";
    cloneBiome->instances.clear();

    Voronoi diagram(this->instances.size(), newArea);
    std::vector<BSpline> subareas = diagram.solve();

//    std::cout << this->instances.size() << " (nb subareas = " << subareas.size() << ")" << std::endl;
    for (size_t i = 0; i < this->instances.size() && i < subareas.size(); i++) {
        std::shared_ptr<BiomeInstance> newChild = this->instances[i]->clone(subareas[i], diagram.pointset[i]);
        newChild->parent = cloneBiome;
        cloneBiome->instances.push_back(newChild);
    }
    return cloneBiome;
}

std::shared_ptr<BiomeInstance> BiomeInstance::getPointInstance(int index)
{
    // Look for each point in the children, find the n-th point by decreasing the index
    for (auto& child : instances) {
        // Uncomment next line to consider any child biome as a 3D point
        // if (child->classname == "point")
        if (child->area)
            index--;
        if (index < 0)
            return child;
    }
    return nullptr;
}

std::string BiomeInstance::getTextureName()
{
    if (this->textureClass.empty() && this->parent != nullptr)
        return this->parent->getTextureName();
    return this->textureClass;
}

std::string BiomeInstance::getInstanceName() const
{
    return this->classname + "#" + std::to_string(this->instanceID);
}

std::vector<std::shared_ptr<BiomeInstance> > BiomeInstance::getPathToChild(std::shared_ptr<BiomeInstance> childToFind)
{
    return childToFind->getPathToParent(this->shared_from_this());
}

std::vector<std::shared_ptr<BiomeInstance> > BiomeInstance::getPathToParent(std::shared_ptr<BiomeInstance> parentToFind)
{
    if (parentToFind == this->shared_from_this()) {
        // If we're looking for me, return just me as a path
        return std::vector<std::shared_ptr<BiomeInstance>>({ this->shared_from_this() });
    } else if (this->parent == nullptr) {
        // If I don't have any parent, we missed the target, return an empty vector, that will be the signal
        return std::vector<std::shared_ptr<BiomeInstance>>();
    } else {
        // Otherwise, ask the parent to look for the target
        std::vector<std::shared_ptr<BiomeInstance>> recursivePath = this->parent->getPathToParent(parentToFind);
        // If we recieved an empty list, the target cannot be find
        if (recursivePath.empty())
            return recursivePath;
        recursivePath.push_back(this->shared_from_this());
        return recursivePath;
    }
}

std::vector<std::shared_ptr<BiomeInstance> > BiomeInstance::getPathTo(std::shared_ptr<BiomeInstance> instanceToFind)
{
    std::vector<std::shared_ptr<BiomeInstance> > path = this->getPathToParent(instanceToFind);
    if (!path.empty()) return path;
    else return this->getPathToChild(instanceToFind);
}

std::vector<std::shared_ptr<BiomeInstance> > BiomeInstance::getAllChildrenBreadthFirst()
{
    std::vector<std::shared_ptr<BiomeInstance> > childrenResults;
    childrenResults.push_back(this->shared_from_this());
    for (auto& child : instances) {
        std::vector<std::shared_ptr<BiomeInstance> > childRes = child->getAllChildrenBreadthFirst();
        childrenResults.insert(childrenResults.end(), childRes.begin(), childRes.end());
    }
    return childrenResults;
}

void BiomeInstance::addInstance(std::shared_ptr<BiomeInstance> newInstance)
{
    BiomeInstance::registerBiomeInstance(newInstance);
    bool replaceWithBlankRegion = false;
    for (auto& instance : instances) {
        if (instance->classname == "point") {
            replaceWithBlankRegion = true;
            Vector3 previousPos = instance->position;
            ShapeCurve previousArea = instance->area;
            instance.swap(newInstance);
            instance->parent = this->shared_from_this();
            instance->position = previousPos;
            instance->area = previousArea;
            instance->updateSubInstances();
            /*break;
            auto previousParent = instance;
//            int previousID = instance->instanceID;

            // Swap the content of the pointers
            *instance = *(newInstance->clone(previousArea));
            // Change the pointer that defines the old parent, so they don't use the old memory address
            for (auto& child : instance->instances) {
                child->parent = instance;
            }

            instance->parent = previousParent;
//            instance->instanceID = previousID;

            instance->completeIfNeeded();*/
            break;
        }
    }
    if (!replaceWithBlankRegion) {
        // No empty biome has been found, we need to create one
        this->instances.push_back(newInstance);
        newInstance->parent = this->shared_from_this();
        Voronoi diagram(this->instances.size(), this->area);
        std::vector<BSpline> subareas = diagram.solve();
        for (size_t i = 0; i < this->instances.size() && i < subareas.size(); i++) {
            this->instances[i]->area = subareas[i];
            this->instances[i]->position = diagram.pointset[i];
            this->instances[i]->updateSubInstances();
        }
    }
}

void BiomeInstance::updateSubInstances()
{
    std::vector<Vector3> newPositions;
    for (auto& instance : instances) {
        if (this->area.contains(instance->position)) {
            newPositions.push_back(instance->position);
        } else {
            std::vector<Vector3> newPossiblePos = this->area.randomPointsInside();
            if (!newPossiblePos.empty())
                newPositions.push_back(newPossiblePos.front());
            else
                newPositions.push_back(Vector3(-1000, -1000));
        }
    }
    std::vector<BSpline> subareas(this->instances.size());
    Voronoi diagram(newPositions, this->area);
    if (this->area.computeArea() > 0) {
        subareas = diagram.solve();
    }
    for (size_t i = 0; i < this->instances.size() && i < subareas.size(); i++) {
        this->instances[i]->area = subareas[i];
        this->instances[i]->position = diagram.pointset[i];
        this->instances[i]->updateSubInstances();
    }

}

std::shared_ptr<BiomeModel> BiomeInstance::toBiomeModel()
{
    // Get the original model used
    std::shared_ptr<BiomeModel> model = this->model->clone();
    // Modify the parameters to make it unique
    model->modelName = this->classname;
    model->textureClass = this->textureClass;

    // Update all the sub-models
    model->modelChildren.clear();
    for (auto& child : this->instances) {
        // Don't consider the empty regions
        if (child->model == nullptr) continue;

        BiomeModelChild modelChild;
        // For now, just consider that all childen are unique...
        modelChild.priorityOffset = child->priorityOffset;
        modelChild.probaAppearence = Probability(1.f, 0.f);
        modelChild.probaQuantity = Probability(1.f, 0.f);
        modelChild.probaSize = Probability(child->area.computeArea(), 0.f);
        // Recursive call
        std::shared_ptr<BiomeModel> childModel = child->toBiomeModel();
        modelChild.model = childModel;
        model->modelChildren.push_back(modelChild);
    }

    return model;
}

void BiomeInstance::registerBiomeInstance(std::shared_ptr<BiomeInstance> biome)
{
    // Don't register simple points
    if (biome->classname == "point") return;
    if (!BiomeInstance::instancedBiomes.empty())
        biome->instanceID = BiomeInstance::instancedBiomes.rbegin()->first + 1;
    else
        biome->instanceID = 0;
    BiomeInstance::instancedBiomes[biome->instanceID] = biome;
}

void BiomeInstance::completeArch()
{
    // We need at least 2 points
    int neededPointsAmount = 2;
    int numberOfGeneratedPoints = neededPointsAmount - getNumberOfPoints();

    Voronoi diagram(numberOfGeneratedPoints + this->instances.size(), this->area);
    std::vector<BSpline> areas = diagram.solve();
    size_t i = 0;
    // Re-place older children
    for (i = 0; i < this->instances.size(); i++) {
        this->instances[i]->area = areas[i];
        this->instances[i]->position = diagram.pointset[i];
    }
    // Add new points
    for (; i < areas.size(); i++) {
        auto child = std::make_shared<BiomeInstance>(BiomeInstance::fromClass("point"));
        instances.push_back(child);
        child->parent = this->shared_from_this();
        child->area = areas[i];
        child->position = diagram.pointset[i];
    }
}

void BiomeInstance::completeTrench()
{
    // We need at least 2 points
    int neededPointsAmount = 2;
    int numberOfGeneratedPoints = neededPointsAmount - getNumberOfPoints();

    Voronoi diagram(numberOfGeneratedPoints + this->instances.size(), this->area);
    std::vector<BSpline> areas = diagram.solve();
    size_t i = 0;
    // Re-place older children
    for (i = 0; i < this->instances.size(); i++) {
        this->instances[i]->area = areas[i];
        this->instances[i]->position = diagram.pointset[i];
    }
    // Add new points
    for (; i < areas.size(); i++) {
        auto child = std::make_shared<BiomeInstance>(BiomeInstance::fromClass("point"));
        instances.push_back(child);
        child->parent = this->shared_from_this();
        child->area = areas[i];
        child->position = diagram.pointset[i];
    }
}

void BiomeInstance::completeArea()
{
    // We need at least 2 points
    int neededPointsAmount = 2;
    int numberOfGeneratedPoints = neededPointsAmount - getNumberOfPoints();

    Voronoi diagram(numberOfGeneratedPoints + this->instances.size(), this->area);
    std::vector<BSpline> areas = diagram.solve();
    size_t i = 0;
    // Re-place older children
    for (i = 0; i < this->instances.size(); i++) {
        this->instances[i]->area = areas[i];
        this->instances[i]->position = diagram.pointset[i];
    }
    // Add new points
    for (; i < areas.size(); i++) {
        auto child = std::make_shared<BiomeInstance>(BiomeInstance::fromClass("point"));
        instances.push_back(child);
        child->parent = this->shared_from_this();
        child->area = areas[i];
        child->position = diagram.pointset[i];
    }
}

int BiomeInstance::getNumberOfPoints()
{
    int nb_points = 0;
    for (auto& child : instances)
        // Uncomment to get only the number of "points" biomes
//        nb_points += 1; //(child->classname == "point" ? 1 : 0);
        if (child->area)
            nb_points++;
    return nb_points;
}
