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

std::shared_ptr<BiomeInstance> BiomeInstance::clone(ShapeCurve newArea)
{
    std::shared_ptr<BiomeInstance> cloneBiome = std::make_shared<BiomeInstance>(*this);
    cloneBiome->area = newArea;
    std::cout << "Previous area : " << this->area << " - new area : " << cloneBiome->area << std::endl;
    std::cout << "Previous instance size : " << this->instances.size() << " - new size : ";
    cloneBiome->instances.clear();

    Voronoi diagram(this->instances.size(), newArea);
    std::vector<BSpline> subareas = diagram.solve();

    std::cout << this->instances.size() << " (nb subareas = " << subareas.size() << ")" << std::endl;
    for (size_t i = 0; i < this->instances.size() && i < subareas.size(); i++) {
        std::shared_ptr<BiomeInstance> newChild = this->instances[i]->clone(subareas[i]);
        newChild->parent = cloneBiome;
        cloneBiome->instances.push_back(newChild);
    }
    return cloneBiome;
}

std::string BiomeInstance::getTextureName()
{
    if (this->textureClass.empty() && this->parent != nullptr)
        return this->parent->getTextureName();
    return this->textureClass;
}

void BiomeInstance::registerBiome(std::shared_ptr<BiomeInstance> biome)
{
    BiomeInstance::instancedBiomes[biome->instanceID] = biome;
}

void BiomeInstance::completeArch()
{
    // We need at least 2 points
    int neededPointsAmount = 2;
    std::vector<Vector3> randomPoints = this->area.randomPointsInside(neededPointsAmount);
    for (size_t i = getNumberOfPoints(); i < randomPoints.size(); i++) {
        auto child = std::make_shared<BiomeInstance>(BiomeInstance::fromClass("point"));
        instances.push_back(child);
        child->position = randomPoints[i];
    }
    if (this->getNumberOfPoints() < neededPointsAmount)
        this->valid = false;
}

void BiomeInstance::completeTrench()
{
    // We need at least 2 points
    int neededPointsAmount = 2;
    std::vector<Vector3> randomPoints = this->area.randomPointsInside(neededPointsAmount);
    for (size_t i = getNumberOfPoints(); i < randomPoints.size(); i++) {
        auto child = std::make_shared<BiomeInstance>(BiomeInstance::fromClass("point"));
        instances.push_back(child);
        child->position = randomPoints[i];
    }
    if (this->getNumberOfPoints() < neededPointsAmount)
        this->valid = false;
}

void BiomeInstance::completeArea()
{
    // We want 10 points, just to get the Hull shape of them
    int neededPointsAmount = 10;
    std::vector<Vector3> randomPoints = this->area.randomPointsInside(neededPointsAmount);
    for (size_t i = getNumberOfPoints(); i < randomPoints.size(); i++) {
        auto child = std::make_shared<BiomeInstance>(BiomeInstance::fromClass("point"));
        instances.push_back(child);
        child->position = randomPoints[i];
    }
    if (this->getNumberOfPoints() < neededPointsAmount)
        this->valid = false;
}

int BiomeInstance::getNumberOfPoints()
{
    int nb_points = 0;
    for (auto& child : instances)
        nb_points += (child->classname == "point" ? 1 : 0);
    return nb_points;
}
