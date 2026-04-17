#ifndef RANDOM_H
#define RANDOM_H

#include <random>
#include "Utils/FastNoiseLit.h"

class random_gen {
public:
    static std::default_random_engine random_generator;
    static FastNoiseLite perlinNoise;
    static float generate(float min, float max) {
        std::uniform_real_distribution<float> distribution(min, max);
        return distribution(random_gen::random_generator);
    }
    static float generate(float max) {
        std::uniform_real_distribution<float> distribution(0.f, max);
        return distribution(random_gen::random_generator);
    }
    static float generate() {
        std::uniform_real_distribution<float> distribution(0.f, 1.f);
        return distribution(random_gen::random_generator);
    }

    static float generate_normal(float mu, float sigma) {
        std::normal_distribution<float> distribution(mu, sigma);
        return distribution(random_gen::random_generator);
    }
    static float generate_normal() {
        std::uniform_real_distribution<float> distribution(0.f, 1.f);
        return distribution(random_gen::random_generator);
    }

    static float generate_perlin(float x, float y = 0, float z = 0) {
        return perlinNoise.GetNoise(x, y, z);
    }

    static int random_int(int min, int max) {
        std::uniform_int_distribution<int> distribution(min, max);
        return distribution(random_gen::random_generator);
    }

    static float generate_fbm(float x, float y = 0, float z = 0, int octaves = 8, float gain = .5f, float lacunarity = 2.f) {
        float value = 0;
        float amplitude = 1.f;
        float frequency = 1.;
        float divisor = 0;
        for (int i = 0; i < octaves; i++) {
            value += random_gen::generate_perlin(x * frequency, y * frequency, z * frequency) * amplitude;
            divisor += amplitude;
            frequency *= lacunarity;
            amplitude *= gain;
        }
        // float divisor = (gain != 0 ? )
        return value / divisor;
    }

    template <typename T>
    static const T& random_choice(const std::vector<T>& vec)
    {
        if (vec.empty()) {
            throw std::runtime_error("random_choice called on empty vector");
        }
        return vec[random_int(0, static_cast<int>(vec.size()) - 1)];
    }
};

#endif // RANDOM_H
