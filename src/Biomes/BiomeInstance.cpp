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

int BiomeInstance::getLevel()
{
    if (!this->parent) {
        return 0;
    } else {
        return this->parent->getLevel() + 1;
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
    cloneBiome->area = newArea;
    if (!newPosition.isValid())
        newPosition = newArea.center();
    cloneBiome->position = newPosition;
    std::cout << "Previous area : " << this->area << " - new area : " << cloneBiome->area << std::endl;
    std::cout << "Previous instance size : " << this->instances.size() << " - new size : ";
    cloneBiome->instances.clear();

    Voronoi diagram(this->instances.size(), newArea);
    std::vector<BSpline> subareas = diagram.solve();

    std::cout << this->instances.size() << " (nb subareas = " << subareas.size() << ")" << std::endl;
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
        if (child->classname == "point")
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

void BiomeInstance::registerBiomeInstance(std::shared_ptr<BiomeInstance> biome)
{
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
        nb_points += (child->classname == "point" ? 1 : 0);
    return nb_points;
}
