#ifndef GLOBALS_H
#define GLOBALS_H

#include <random>
class random_gen {
public:
    static std::default_random_engine random_generator;
    static float generate(float min = 0.0, float max = 1.0) {
        std::uniform_real_distribution<float> distribution(min, max);
        return distribution(random_gen::random_generator);
    }
};
#endif // GLOBALS_H
