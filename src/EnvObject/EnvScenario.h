#ifndef ENVSCENARIO_H
#define ENVSCENARIO_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "Utils/Globals.h"

#include "DataStructure/Matrix3.h"


class Scenario;

struct ScenariosObject {
    ScenariosObject(std::string objectName, float proba, int amount = -1) : objectName(objectName), probabilityPerStep(proba), amountRequired(amount)
    {}
    std::string objectName;
    float probabilityPerStep;
    int amountRequired;

    float normalizedProba;
};

class WaterLevelEvent {
public:
    WaterLevelEvent(float amount, float startTime, float endTime)
        : amount(amount), startTime(startTime), endTime(endTime)
    {}

    float applyAt(const Vector3& pos, float progress);
    float apply(float progress, const Vector3& dimensions = Vector3());

    float amount;
    float startTime;
    float endTime;
};

class SubsidenceEvent {
public:
    SubsidenceEvent(const Vector3& position, float amount, float sigma, float startTime, float endTime)
        : position(position), amount(amount), sigma(sigma), startTime(startTime), endTime(endTime)
    {}

    float applyAt(const Vector3& pos, float progress);
    GridF apply(float progress, const Vector3& dimensions);

    Vector3 position;
    float amount;
    float sigma;
    float startTime;
    float endTime;
};

class StormEvent {
public:
    StormEvent(const Vector3& pos, const Vector3& dir, float sigma, float startTime, float endTime)
        : position(pos), direction(dir), sigma(sigma), startTime(startTime), endTime(endTime)
    {}

    Vector3 applyAt(const Vector3& pos, float progress);
    GridV3 apply(float progress, const Vector3 &dimensions);

    Vector3 position;
    Vector3 direction;
    float sigma;
    float startTime;
    float endTime;
};

class TectonicEvent {
public:
    TectonicEvent(Vector3 direction, float sigma, float startTime, float endTime);

    float applyAt(const Vector3& pos, float progress);
    GridF apply(float progress, const Vector3& dimensions);

    Vector3 direction;
    float sigma;
    float startTime;
    float endTime;

    GridF initialState;
};

class WarmingEvent {
public:
    WarmingEvent(float amount, float startTime, float endTime)
        : amount(amount), startTime(startTime), endTime(endTime)
    {}

    float applyAt(const Vector3& pos, float progress);
    float apply(float progress, const Vector3& dimensions = Vector3());

    float amount;
    float startTime;
    float endTime;

    GridF initialState;
};

class Scenario {
public:
    void addObject(std::string name, float proba, int amount = -1);

    std::vector<ScenariosObject> objects;

    std::vector<WaterLevelEvent> waterLevelEvents;
    std::vector<SubsidenceEvent> subsidenceEvents;
    std::vector<StormEvent> stormEvents;
    std::vector<TectonicEvent> tectonicEvents;
    std::vector<WarmingEvent> warmingEvents;

    ScenariosObject nextObject();
    std::vector<ScenariosObject> nextObjects();

    float computeWaterLevel(const Vector3 &dimensions = Vector3());
    GridF computeSubsidence(const Vector3 &dimensions);
    GridV3 computeStorm(const Vector3 &dimensions);
    GridF computeTectonic(const Vector3 &dimensions);
    float computeWarming(const Vector3 &dimensions = Vector3());

    float currentTime() const;

    float duration;
    float dt = 1.f;
    float startTime = 0;

    float waterLevel = 0;
    float subsidenceLevel = 1.f;

    GridF subsidence;
    GridV3 storms;
    GridF tectonic;

    bool finished() const;
};


#endif // ENVSCENARIO_H
