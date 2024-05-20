#include "EnvScenario.h"

#include "EnvObject/EnvObject.h"

void Scenario::addObject(std::string name, float proba, int amount)
{
    for (int i = objects.size() - 1; i >= 0; i--) {
        if (objects[i].objectName == name)
            objects.erase(objects.begin() + i);
    }
    objects.push_back(ScenariosObject(name, proba, amount));

    // Recompute the normalized probability
    float sumProbas = 0;
    for (auto& object : objects) {
        sumProbas += object.probabilityPerStep;
    }
    for (auto& object : objects) {
        object.normalizedProba = object.probabilityPerStep / sumProbas;
    }
}

ScenariosObject Scenario::nextObject()
{
    float sumProbas = 0;
    std::map<size_t, float> probas;
    for (size_t i = 0; i < objects.size(); i++) {
        auto& requiredObjects = objects[i];
        int objCount = 0;
        for (auto& obj : EnvObject::instantiatedObjects) {
            if (obj->name == requiredObjects.objectName)
                objCount ++;
        }
        if (objCount < requiredObjects.amountRequired || requiredObjects.amountRequired < 0) {
            sumProbas += requiredObjects.normalizedProba;
            probas[i] = requiredObjects.normalizedProba;
        }
    }



    float isNextObject = random_gen::generate(sumProbas);
    for (auto& [i, proba] : probas) {
        // std::cout << objects[i].objectName << " -> " << isNextObject << " .. " << std::flush;
        isNextObject -= proba;

        if (isNextObject <= 0) return objects[i];
    }

    std::cerr << "WTF, should not be here..." << std::endl;
    return objects[int(random_gen::generate(objects.size()))];
}

bool Scenario::finished() const
{
    if (duration > 0 && EnvObject::currentTime >= startTime + duration)
        return true;

    for (auto& requiredObjects : objects) {
        if (requiredObjects.amountRequired < 0) continue;
        int objCount = 0;
        for (auto& obj : EnvObject::instantiatedObjects) {
            if (obj->name == requiredObjects.objectName)
                objCount ++;
        }
        if (objCount < requiredObjects.amountRequired) {
            return false;
        }
    }
    return true;
}


ScenarioEvent::ScenarioEvent(std::string typeName, float amount, float startTime, float endTime)
    : typeName(toLower(typeName)), amount(amount), startTime(startTime), endTime(endTime)
{
    if (typeName == "waterlevel") {
        type = WATER_LEVEL;
    } else if (typeName == "material_deposition") {
        type = MATERIAL_DEPOSITION;
    }
}

