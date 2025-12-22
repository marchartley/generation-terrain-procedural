#include "PSO.h"

PSO::PSO() {}

std::pair<Vector3, float> PSO::findHighest(size_t nbParticles, const AABBox &domain, std::function<float (const Vector3 &)> func, int maxIterations)
{
    std::vector<PSOParticle> particles(nbParticles);
    float currentBestValue = std::numeric_limits<float>::lowest();
    size_t bestParticle = 0;
    // Init
    for (size_t i = 0; i < nbParticles; i++) {
        particles[i].pos = Vector3::random(domain.mini, domain.maxi);
        particles[i].currentBestPos = particles[i].pos;
        float eval = func(particles[i].pos);
        particles[i].currentBestVal = eval;
        if (eval > currentBestValue) {
            currentBestValue = eval;
            bestParticle = i;
        }
    }

    Vector3 currentBestPos = particles[bestParticle].pos;

    float a = 0.8; // damping weight
    float b = 0.8; // cognitive weight
    float b_hat = 1; // social weight
    float c = 0.1;

    for (int iter = 0; iter < maxIterations; iter++) {
        #pragma omp parallel for
        for (size_t i = 0; i < nbParticles; i++) {
            auto& particle = particles[i];
            float r = random_gen::generate(); // randomness
            float r_hat = random_gen::generate(); // randomness

            // Cognitive component: Tendency towards personal best
            Vector3 cognitiveVelocity = (particle.currentBestPos - particle.pos) * (b * r);

            // Social component: Tendency towards global best
            Vector3 socialVelocity = (currentBestPos - particle.pos) * (b_hat * r_hat);

            // Update velocity with inertia
            particle.vel = (particle.vel * a) + (cognitiveVelocity + socialVelocity);

            // Update position
            particle.pos += particle.vel * c;

            // Update the evaluation based on the new position
            float currentEval = func(particle.pos);
            if (currentEval > particle.currentBestVal) {
                particle.currentBestVal = currentEval;
                particle.currentBestPos = particle.pos;
            }
        }
        for (size_t i = 0; i < nbParticles; i++) {
            auto& particle = particles[i];
            if (particle.currentBestVal > currentBestValue) {
                currentBestValue = particle.currentBestVal;
                currentBestPos = particle.pos;
            }
        }
    }
    return {currentBestPos, currentBestValue};
}

std::vector<std::vector<std::pair<Vector3, float>>> PSO::findHighestAndTrack(size_t nbParticles, const AABBox &domain, std::function<float (const Vector3 &)> func, int maxIterations)
{
    std::vector<PSOParticle> particles(nbParticles);
    std::vector<std::vector<std::pair<Vector3, float>>> trajectories(nbParticles, std::vector<std::pair<Vector3, float>>(maxIterations));
    float currentBestValue = std::numeric_limits<float>::lowest();
    size_t bestParticle = 0;
    // Init
    for (size_t i = 0; i < nbParticles; i++) {
        particles[i].pos = Vector3::random(domain.mini, domain.maxi);
        particles[i].currentBestPos = particles[i].pos;
        float eval = func(particles[i].pos);
        particles[i].currentBestVal = eval;
        if (eval > currentBestValue) {
            currentBestValue = eval;
            bestParticle = i;
        }
    }

    Vector3 currentBestPos = particles[bestParticle].pos;

    float a = 0.8; // damping weight
    float b = 0.8; // cognitive weight
    float b_hat = 1; // social weight
    float c = 0.1;

    for (int iter = 0; iter < maxIterations; iter++) {
        #pragma omp parallel for
        for (size_t i = 0; i < nbParticles; i++) {
            auto& particle = particles[i];
            float r = random_gen::generate(); // randomness
            float r_hat = random_gen::generate(); // randomness

            // Cognitive component: Tendency towards personal best
            Vector3 cognitiveVelocity = (particle.currentBestPos - particle.pos) * (b * r);

            // Social component: Tendency towards global best
            Vector3 socialVelocity = (currentBestPos - particle.pos) * (b_hat * r_hat);

            // Update velocity with inertia
            particle.vel = (particle.vel * a) + (cognitiveVelocity + socialVelocity);

            // Update position
            particle.pos += particle.vel * c;

            // Update the evaluation based on the new position
            float currentEval = func(particle.pos);
            if (currentEval > particle.currentBestVal) {
                particle.currentBestVal = currentEval;
                particle.currentBestPos = particle.pos;
            }
            trajectories[i][iter] = {particle.pos, currentEval};
        }
        for (size_t i = 0; i < nbParticles; i++) {
            auto& particle = particles[i];
            if (particle.currentBestVal > currentBestValue) {
                currentBestValue = particle.currentBestVal;
                currentBestPos = particle.pos;
            }
        }
    }
    return trajectories;
}
