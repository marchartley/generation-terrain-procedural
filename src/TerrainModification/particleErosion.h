#ifndef PARTICLEEROSION_H
#define PARTICLEEROSION_H

#include "DataStructure/Vector3.h"
#include "DataStructure/Matrix3.h"


struct ErosionParticle {
    Vector3 pos;
    Vector3 dir;
    Vector3 velocity;
    float capacity;
    float maxCapacity;
    Vector3 force;
    float density;
    float volume;
    float radius;
    float mass;
};

std::vector<std::vector<Vector3>> erosion();
void _initializeParticle(ErosionParticle& particle, Vector3& position, Vector3& velocity, float radius, float density, float initialCapacity, float maxCapacity);

#endif // PARTICLEEROSION_H
