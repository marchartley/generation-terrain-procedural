#ifndef PARTICLE_H
#define PARTICLE_H

#include "DataStructure/Vector3.h"

class Particle {
public:
    Vector3 position;
    Vector3 velocity;
    Vector3 force;
    float density = 0;
    float pressure = 0;
    float mass = 0;
    float smoothingRadius = 0;
    float gasConstant = 0;
    float restDensity = 0;
    float viscosity = 0;
    int isGhost = 0;

    int index = -1;
};

#endif // PARTICLE_H
