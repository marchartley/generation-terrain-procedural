#ifndef NEAT_H
#define NEAT_H
#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Utils/Random.h"

enum class NodeKind { Input, Hidden, Output };
enum class AggregationKind { Sum, Product };
enum class ActivationFunction { Sigmoid, ReLU, Tanh, Linear, Cos, Sqr, Sqrt };


// ----------------------------
// Activation
// ----------------------------
float sigmoid(float x);

float reLU(float x);

// ----------------------------
// Global node ID manager
// ----------------------------
class NodeIDGenerator {
public:
    explicit NodeIDGenerator(int start = 0);

    int new_id();

    NodeIDGenerator copy() const;

private:
    int next_id;
};

// ----------------------------
// Gene definitions
// ----------------------------
struct NodeGene {
    NodeGene() {}
    NodeGene(int id, NodeKind kind, AggregationKind aggregation, ActivationFunction activation, float bias);
    int id;
    NodeKind kind;       // "input", "hidden", "output"
    AggregationKind aggregation;
    ActivationFunction activation;
    float bias = 0.0;
};

struct ConnectionGene {
    int in_node;
    int out_node;
    float weight;
    bool enabled;
};

// ----------------------------
// Genome
// ----------------------------
class Genome {
public:
    std::map<int, NodeGene> nodes;
    std::vector<ConnectionGene> connections;
    float fitness = -1000000.0;
    int generation = 0;
    std::vector<int> input_node_ids;
    std::vector<int> output_node_ids;
    NodeIDGenerator node_id_gen;

    Genome(int nb_inputs = 0, int nb_outputs = 0);

    inline static std::vector<AggregationKind> available_aggregation_funcs();
    inline static std::vector<ActivationFunction> available_activation_funcs();

    static float apply_activation(ActivationFunction func, float x);
    static std::string get_activation_name(ActivationFunction func);

    void get_input_output_node_ids();

    Genome copy() const;

    // ----------------------------
    // Network execution
    // ----------------------------
    std::vector<int> topological_order();

    std::vector<int> get_useful_nodes() const;

    std::vector<float> activate_network(const std::vector<float>& input_values);

    // ----------------------------
    // Mutation
    // ----------------------------
    void mutate_weights(
        float weight_mutation_rate,
        float max_weight_perturb,
        float activation_function_rate
        );

    bool connection_exists(int in_node, int out_node) const;

    bool would_create_cycle(int in_node, int out_node) const;

    void add_connection_mutation(float adding_rate);

    void remove_connection_mutation(float removal_rate);

    void add_node_mutation();

    void remove_node_mutation();

    void mutate(
        float weight_mutation_rate,
        float add_connection_rate,
        float add_node_rate,
        float remove_connection_rate,
        float remove_node_rate,
        float activation_function_rate,
        float max_weight_perturb
        );

    std::string get_symbolic_formula(
        const std::vector<std::string>& input_names,
        std::optional<std::vector<int>> output_ids = std::nullopt
        ) const;

    std::string summarize() const;

protected:
    std::vector<int> cached_topological_order;
    bool dirty_topological_order = true;
};

// ----------------------------
// Population
// ----------------------------
class Population {
public:
    int population_size = 50;
    int generations = 10;

    float elite_fraction = 0.1;
    float weight_mutation_rate = 0.8;
    float add_connection_rate = 0.5;
    float add_node_rate = 0.5;
    float remove_connection_rate = 0.1;
    float remove_node_rate = 0.3;
    float activation_function_rate = 0.2;
    float network_size_loss_factor = 1.01;

    float max_weight_perturb = 0.5;

    int current_generation = 0;

    std::vector<Genome> genomes;
    Genome current_best_genome;

    Population& init(int nb_inputs, int nb_outputs);

    std::vector<Genome> reproduce();

    template <typename FitnessFunction>
    Genome run_once(FitnessFunction fitness_function) {
        current_generation += 1;
        #pragma omp parallel for
        for (auto& genome : genomes) {
            genome.fitness = fitness_function(genome);
        }

        std::sort(genomes.begin(), genomes.end(),
                  [](const Genome& a, const Genome& b) {
                      return a.fitness > b.fitness;
                  });

        const Genome& best = genomes[0];

        if (best.fitness > current_best_genome.fitness) {
            current_best_genome = best.copy();
        }

        reproduce();
        return current_best_genome;
    }
    template <typename FitnessFunction>
    std::vector<Genome> run(FitnessFunction fitness_function) {
        std::vector<Genome> history;

        if (genomes.empty()) {
            throw std::runtime_error("Population is empty. Call init() first.");
        }

        for (int generation = 0; generation < generations; ++generation) {
            history.push_back(run_once(fitness_function));
        }

        return history;
    }
};

// ----------------------------
// Example usage
// ----------------------------
// int main() {
//     Population pop;
//     pop.population_size = 20;
//     pop.generations = 5;
//     pop.init(2, 1);

//     auto fitness_function = [](const Genome& genome) -> float {
//         // Dummy example fitness
//         std::vector<float> out = genome.activate_network({1.0, 2.0});
//         return out.empty() ? 0.0 : out[0];
//     };

//     std::vector<Genome> history = pop.run(fitness_function);

//     for (const auto& best : history) {
//         std::cout << best.summarize() << "\n";
//     }

//     if (!history.empty()) {
//         std::cout << "Formula: "
//                   << history.back().get_symbolic_formula({"x1", "x2"})
//                   << "\n";
//     }

//     return 0;
// }
#endif // NEAT_H
