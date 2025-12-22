#ifndef PSO_H
#define PSO_H

#include "DataStructure/Vector3.h"
#include <utility>
#include <functional>
#include <cstddef>

class PSOParticle;

class PSO
{
public:
    PSO();
    static std::pair<Vector3, float> findHighest(size_t nbParticles, const AABBox& domain, std::function<float(const Vector3&)> func, int maxIterations = 100);
    static std::vector<std::vector<std::pair<Vector3, float> > > findHighestAndTrack(size_t nbParticles, const AABBox& domain, std::function<float(const Vector3&)> func, int maxIterations = 100);
};

class PSOParticle {
public:
    Vector3 pos;
    Vector3 vel;
    Vector3 currentBestPos;
    float currentBestVal;
};
#endif // PSO_H
