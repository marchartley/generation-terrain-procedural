#ifndef ENVSCENARIO_H
#define ENVSCENARIO_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "Utils/Globals.h"

struct ScenariosObject {
    ScenariosObject(std::string objectName, float proba, int amount = -1) : objectName(objectName), probabilityPerStep(proba), amountRequired(amount)
    {}
    std::string objectName;
    float probabilityPerStep;
    int amountRequired;

    float normalizedProba;
};

struct ScenarioEvent {
    enum Type {
        WATER_LEVEL,
        SUBSIDENCE,
        MATERIAL_DEPOSITION
    };

    ScenarioEvent(std::string typeName, float amount, float startTime, float endTime);

    Type type;
    std::string typeName;
    float amount;
    float startTime;
    float endTime;
};

class Scenario {
public:
    void addObject(std::string name, float proba, int amount = -1);
    std::vector<ScenariosObject> objects;
    ScenariosObject nextObject();


    float duration;
    float dt = 1.f;
    float startTime = 0;

    bool finished() const;
};


#endif // ENVSCENARIO_H
