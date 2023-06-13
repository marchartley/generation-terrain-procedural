#include "SPHSimulation.h"

float SPHSimulation::distance(Particle &p1, Particle &p2) {
    return (Vector3(p1.x, p1.y, p1.z) - Vector3(p2.x, p2.y, p2.z)).norm();
}
