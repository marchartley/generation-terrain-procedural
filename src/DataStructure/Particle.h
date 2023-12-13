#ifndef PARTICLE_H
#define PARTICLE_H

#include "DataStructure/Vector3.h"

class Particle {
public:
    Vector3 position;
    Vector3 velocity;
    Vector3 force;
    float density;
    float pressure;
    float mass;
    float smoothingRadius;
    float gasConstant;
    float restDensity;
    float viscosity;
    int isGhost;

    int index;
};

#endif // PARTICLE_H
