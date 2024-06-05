#include "EnvScenario.h"

#include "EnvObject/EnvObject.h"

#include "Utils/Delaunay.h"

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
        isNextObject -= proba;

        if (isNextObject <= 0) return objects[i];
    }

    std::cerr << "WTF, should not be here..." << std::endl;
    return objects[int(random_gen::generate(objects.size()))];
}

std::vector<ScenariosObject> Scenario::nextObjects()
{
    std::vector<ScenariosObject> results;
    for (auto& possibleObj : this->objects) {
        float proba = possibleObj.probabilityPerStep;
        int objCount = 0;
        for (auto& obj : EnvObject::instantiatedObjects) {
            if (obj->name == possibleObj.objectName)
                objCount ++;
        }
        for (int i = 0; i < dt && (objCount < possibleObj.amountRequired || possibleObj.amountRequired < 0); i++) {
            // std::cout << "Generate " << possibleObj.objectName << "? " << possibleObj.amountRequired << " / " << objCount << std::endl;
            if (random_gen::generate() < proba) {
                results.push_back(possibleObj);
                objCount ++;
            }
        }
    }
    return results;
}

float Scenario::computeWaterLevel(const Vector3 &dimensions)
{
    float time = currentTime();
    float water = this->waterLevel;
    for (auto& event : waterLevelEvents) {
        // Should be 1 for all finished events
        float progress = 1;
        if (time < event.endTime)
            progress = interpolation::linear(time, event.startTime, event.endTime);
        water += event.apply(progress);
    }
    return water;
}

GridF Scenario::computeSubsidence(const Vector3& dimensions)
{
    float time = currentTime();
    GridF subsidence(dimensions, this->subsidenceLevel);
    for (auto& event : subsidenceEvents) {
        // Should be 1 for all finished events
        float progress = interpolation::linear(time, event.startTime, event.endTime);
        subsidence += event.apply(progress, dimensions);
    }
    subsidence.iterateParallel([&](size_t i) {
        subsidence[i] = std::clamp(subsidence[i], 0.f, 1.f);
    });
    return subsidence;
}

GridV3 Scenario::computeStorm(const Vector3 &dimensions)
{
    float time = currentTime();
    GridV3 storm(dimensions, Vector3());
    for (auto& event : stormEvents) {
        // Should be 1 for all finished events
        float progress = interpolation::linear(time, event.startTime, event.endTime);
        storm += event.apply(progress, dimensions);
    }
    return storm;
}

GridF Scenario::computeTectonic(const Vector3 &dimensions)
{
    float time = currentTime();
    GridF tecto(dimensions, 0.f);
    for (auto& event : tectonicEvents) {
        // Should be 1 for all finished events
        float progress = interpolation::linear(time, event.startTime, event.endTime);
        tecto += event.apply(progress, dimensions);
    }
    return tecto;
}

float Scenario::computeWarming(const Vector3 &dimensions)
{
    float time = currentTime();
    float warm = 1.f;
    for (auto& event : warmingEvents) {
        // Should be 1 for all finished events
        float progress = interpolation::linear(time, event.startTime, event.endTime);
        warm += event.apply(progress, dimensions);
    }
    return warm;
}

float Scenario::currentTime() const
{
    return EnvObject::currentTime - startTime;
}

bool Scenario::finished() const
{
    float time = currentTime();
    if (duration > 0 && time >= duration)
        return true;

    bool finish = true;
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

    for (auto& event : waterLevelEvents) {
        if (event.endTime > time)
            return false;
    }
    for (auto& event : stormEvents) {
        if (event.endTime > time)
            return false;
    }
    for (auto& event : tectonicEvents) {
        if (event.endTime > time)
            return false;
    }
    for (auto& event : subsidenceEvents) {
        if (event.endTime > time)
            return false;
    }
    return true;
}

float WaterLevelEvent::applyAt(const Vector3 &pos, float progress)
{
    if (progress < 0) return 0;
    return amount * std::clamp(progress, 0.f, 1.f);
}

float WaterLevelEvent::apply(float progress, const Vector3 &dimensions)
{
    return applyAt(Vector3(), progress);
}

float SubsidenceEvent::applyAt(const Vector3 &pos, float progress)
{
    if (progress < 0) return 0;
    if (sigma > 0) {
        return -normalizedGaussian(this->sigma, (pos - position).norm2()) * amount * std::clamp(progress, 0.f, 1.f);
    } else {
        return -amount * std::clamp(progress, 0.f, 1.f);
    }
}

GridF SubsidenceEvent::apply(float progress, const Vector3 &dimensions)
{
    GridF result(dimensions, 0.f);
    if (progress <= 0) return result;

    result.iterateParallel([&](const Vector3& pos) {
        result(pos) = applyAt(pos, progress);
    });
    return result;
}

Vector3 StormEvent::applyAt(const Vector3 &pos, float progress)
{
    if (progress <= 0) return Vector3();
    return direction * normalizedGaussian(this->sigma, (pos - position).norm2()) * std::clamp(progress, 0.f, 1.f) * normalizedGaussian(.3f, (progress - 0.5f)*(progress - 0.5f));
}

GridV3 StormEvent::apply(float progress, const Vector3& dimensions)
{
    GridV3 result(dimensions, Vector3());
    if (progress <= 0) return result;

    result.iterateParallel([&](const Vector3& pos) {
        result(pos) = applyAt(pos, progress);
    });
    return result;

}

TectonicEvent::TectonicEvent(Vector3 direction, float sigma, float startTime, float endTime)
    : direction(direction), sigma(sigma), startTime(startTime), endTime(endTime)
{
    this->initialState = GridF(100, 100, 1);

    Voronoi voro(20, Vector3(-10, -10, 0), Vector3(110, 110, 0));
    voro.solve(false);
    #pragma omp parallel for
    for (auto& area : voro.areas) {
        for (int i = 0; i < area.size(); i++) {
            Vector3 p0 = area[i];
            Vector3 p1 = area[i+1];
            BSpline path = BSpline({p0, p1}).resamplePoints(10);
            Vector3 dir = (p0 - p1);
            dir = Vector3(dir.y, dir.x).normalize();
            for (auto& p : path) {
                p += random_gen::generate_normal() * 2.f * dir;
            }
            for (auto& p : path.getPath(50)) {
                initialState(p) = 1.f;
            }
        }
    }
    /*
    Delaunay delaunay(voro);

    #pragma omp parallel for
    for (auto& n : delaunay.graph.nodes) {
        for (auto& [nn, dist] : n->neighbors) {
            BSpline path = BSpline({n->pos, nn->pos}).resamplePoints(10);
            Vector3 dir = (n->pos - nn->pos);
            dir = Vector3(dir.y, dir.x).normalize();
            for (auto& p : path) {
                p += random_gen::generate_normal() * 2.f * dir;
            }
            for (auto& p : path.getPath(100)) {
                initialState(p) = 1.f;
            }
        }
    }*/
}

float TectonicEvent::applyAt(const Vector3 &pos, float progress)
{
    if (progress <= 0) return 0.f;
    float result = 0.f;

    for (int i = 0; i < direction.length() * 2.f; i++) {
        float t = float(i) / (direction.length() * 2.f);
        result += initialState(pos + direction * t) / (2 * direction.length());
    }
    return result * sigma * std::clamp(progress, 0.f, 1.f);
}

GridF TectonicEvent::apply(float progress, const Vector3 &dimensions)
{
    GridF result(dimensions, 0.f);
    if (progress <= 0) return result;

    result.iterateParallel([&](const Vector3& pos) {
        result(pos) = applyAt(pos, progress);
    });
    return result;
}


float WarmingEvent::applyAt(const Vector3 &pos, float progress)
{
    if (progress <= 0) return 0.f;
    return std::clamp(progress, 0.f, 1.f) * amount;
}

float WarmingEvent::apply(float progress, const Vector3 &dimensions)
{
    return this->applyAt(Vector3(), progress);
}
