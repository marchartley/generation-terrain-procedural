#include "Neat.h"

#include "Utils/Utils.h"

NodeIDGenerator::NodeIDGenerator(int start) : next_id(start) {}

int NodeIDGenerator::new_id() {
    return next_id++;
}

NodeIDGenerator NodeIDGenerator::copy() const {
    return NodeIDGenerator(next_id);
}

Genome::Genome(int nb_inputs, int nb_outputs)
    : node_id_gen(nb_inputs + nb_outputs) {
    for (int node_id = 0; node_id < nb_inputs; ++node_id) {
        nodes[node_id] = NodeGene(node_id, NodeKind::Input, AggregationKind::Sum, ActivationFunction::Linear, 0.0);
    }

    for (int node_id = nb_inputs; node_id < nb_inputs + nb_outputs; ++node_id) {
        nodes[node_id] = NodeGene(node_id, NodeKind::Output, AggregationKind::Sum, ActivationFunction::Linear, 0.0);
    }

    get_input_output_node_ids();

    for (int in_node : input_node_ids) {
        for (int out_node : output_node_ids) {
            connections.push_back(ConnectionGene{
                in_node,
                out_node,
                random_gen::generate(-1.0, 1.0),
                true
            });
        }
    }
}

std::vector<AggregationKind> Genome::available_aggregation_funcs()
{
    return {
        AggregationKind::Sum,
            AggregationKind::Product
    };
}

std::vector<ActivationFunction> Genome::available_activation_funcs() {
    return {
        // ActivationFunction::Cos,
        ActivationFunction::Linear,
        // ActivationFunction::ReLU,
        // ActivationFunction::Sigmoid,
        // ActivationFunction::Tanh,
        // ActivationFunction::Sqr,
        // ActivationFunction::Sqrt
    };
}

float Genome::apply_activation(ActivationFunction func, float x) {
    switch (func) {
    case ActivationFunction::Cos: return std::cos(x);
    case ActivationFunction::Linear: return x;
    case ActivationFunction::ReLU: return reLU(x);
    case ActivationFunction::Sigmoid: return sigmoid(x);
    case ActivationFunction::Tanh: return std::tanh(x);
    case ActivationFunction::Sqr: return x * x;
    case ActivationFunction::Sqrt: return std::sqrt(x);
    }
    throw std::runtime_error("Unknown activation function.");
}

std::string Genome::get_activation_name(ActivationFunction func)
{
    switch (func) {
    case ActivationFunction::Cos: return "cos";
    case ActivationFunction::Linear: return "linear";
    case ActivationFunction::ReLU: return "reLU";
    case ActivationFunction::Sigmoid: return "sigmoid";
    case ActivationFunction::Tanh: return "tanh";
    case ActivationFunction::Sqr: return "pow2";
    case ActivationFunction::Sqrt: return "sqrt";
    }
    throw std::runtime_error("Unknown activation function.");
}

void Genome::get_input_output_node_ids() {
    input_node_ids.clear();
    output_node_ids.clear();

    for (const auto& [id, node] : nodes) {
        if (node.kind == NodeKind::Input) input_node_ids.push_back(id);
        if (node.kind == NodeKind::Output) output_node_ids.push_back(id);
    }
}

Genome Genome::copy() const {
    Genome new_genome;
    new_genome.nodes = nodes;
    new_genome.input_node_ids = input_node_ids;
    new_genome.output_node_ids = output_node_ids;
    new_genome.node_id_gen = node_id_gen.copy();
    new_genome.connections = connections;
    new_genome.fitness = fitness;
    new_genome.generation = generation;
    new_genome.cached_topological_order = cached_topological_order;
    new_genome.dirty_topological_order = dirty_topological_order;
    return new_genome;
}

std::vector<int> Genome::topological_order() {
    if (!dirty_topological_order) return this->cached_topological_order;

    // auto useful_nodes = get_useful_nodes();

    std::map<int, int> incoming_count;
    std::map<int, std::vector<int>> outgoing;

    for (const auto& [node_id, _] : nodes) {
        // if (!isIn(node_id, useful_nodes)) continue;
        incoming_count[node_id] = 0;
        outgoing[node_id] = {};
    }

    for (const auto& conn : connections) {
        if (!conn.enabled) continue;
        outgoing[conn.in_node].push_back(conn.out_node);
        incoming_count[conn.out_node] += 1;
    }

    std::vector<int> queue;
    for (const auto& [node_id, count] : incoming_count) {
        if (count == 0) queue.push_back(node_id);
    }

    std::vector<int> order;
    std::size_t idx = 0;

    while (idx < queue.size()) {
        int node = queue[idx++];
        order.push_back(node);

        for (int neighbor : outgoing[node]) {
            incoming_count[neighbor] -= 1;
            if (incoming_count[neighbor] == 0) {
                queue.push_back(neighbor);
            }
        }
    }

    if (order.size() != nodes.size()) {
        throw std::runtime_error("Cycle detected. This demo only supports acyclic networks.");
    }
    this->cached_topological_order = order;
    this->dirty_topological_order = false;
    return order;
}

std::vector<int> Genome::get_useful_nodes() const {
    std::set<int> used_nodes;
    std::set<int> queue(output_node_ids.begin(), output_node_ids.end());

    while (!queue.empty()) {
        int node = *queue.begin();
        queue.erase(queue.begin());

        used_nodes.insert(node);

        for (const auto& conn : connections) {
            if (!conn.enabled) continue;
            if (conn.out_node == node) {
                queue.insert(conn.in_node);
            }
        }
    }

    return std::vector<int>(used_nodes.begin(), used_nodes.end());
}

std::vector<float> Genome::activate_network(const std::vector<float> &input_values) {
    std::map<int, float> values;
    for (const auto& [node_id, _] : nodes) {
        values[node_id] = 0.0;
    }

    if (input_values.size() != input_node_ids.size()) {
        throw std::runtime_error("Input length does not match number of input nodes.");
    }

    for (std::size_t i = 0; i < input_node_ids.size(); ++i) {
        values[input_node_ids[i]] = input_values[i];
    }

    std::vector<int> order = topological_order();

    std::map<int, std::vector<ConnectionGene>> incoming_map;
    for (const auto& [node_id, _] : nodes) {
        incoming_map[node_id] = {};
    }

    for (const auto& conn : connections) {
        if (conn.enabled) {
            incoming_map[conn.out_node].push_back(conn);
        }
    }

    for (int node_id : order) {
        const NodeGene& node = nodes.at(node_id);
        float total = (node.aggregation == AggregationKind::Sum ? 0.0 : 1.0);

        if (node.kind == NodeKind::Input) {
            total = values[node_id];
        }

        else if (node.aggregation == AggregationKind::Sum) {
            for (const auto& conn : incoming_map[node_id])
                total += values[conn.in_node] * conn.weight;
        } else if (node.aggregation == AggregationKind::Product) {
            for (const auto& conn : incoming_map[node_id])
                total *= (values[conn.in_node] * conn.weight);
        }

        values[node_id] = apply_activation(node.activation, total) + node.bias;
    }

    std::vector<float> outputs;
    for (int node_id : output_node_ids) {
        outputs.push_back(values[node_id]);
    }

    return outputs;
}

void Genome::mutate_weights(float weight_mutation_rate, float max_weight_perturb, float activation_function_rate) {
    for (auto& conn : connections) {
        if (random_gen::generate() < weight_mutation_rate) {
            conn.weight += random_gen::generate(-max_weight_perturb, max_weight_perturb);
        }
    }

    auto acts = available_activation_funcs();

    for (auto& [_, node] : nodes) {
        if (random_gen::generate() < weight_mutation_rate) {
            node.bias += random_gen::generate(-max_weight_perturb, max_weight_perturb);
        }

        if (node.kind == NodeKind::Input || node.kind == NodeKind::Output) continue;

        if (random_gen::generate() < activation_function_rate) {
            node.activation = random_gen::random_choice(acts);
        }
    }
}

bool Genome::connection_exists(int in_node, int out_node) const {
    for (const auto& conn : connections) {
        if (conn.in_node == in_node && conn.out_node == out_node) {
            return true;
        }
    }
    return false;
}

bool Genome::would_create_cycle(int in_node, int out_node) const {
    if (in_node == out_node) {
        return true;
    }

    std::map<int, std::vector<int>> adjacency;
    for (const auto& [node_id, _] : nodes) {
        adjacency[node_id] = {};
    }

    for (const auto& conn : connections) {
        if (conn.enabled) {
            adjacency[conn.in_node].push_back(conn.out_node);
        }
    }

    std::vector<int> stack = {out_node};
    std::set<int> visited;

    while (!stack.empty()) {
        int node = stack.back();
        stack.pop_back();

        if (node == in_node) {
            return true;
        }
        if (visited.count(node)) {
            continue;
        }

        visited.insert(node);
        for (int neighbor : adjacency[node]) {
            stack.push_back(neighbor);
        }
    }

    return false;
}

void Genome::add_connection_mutation(float adding_rate) {
    std::vector<int> node_ids;
    for (const auto& [id, _] : nodes) {
        node_ids.push_back(id);
    }
    std::sort(node_ids.begin(), node_ids.end());

    std::vector<std::pair<int, int>> possible_pairs;

    for (int in_node : node_ids) {
        for (int out_node : node_ids) {
            if (in_node >= out_node) continue;

            const auto& in_kind = nodes.at(in_node).kind;
            const auto& out_kind = nodes.at(out_node).kind;

            if (out_kind == NodeKind::Input) continue;
            if (in_kind == NodeKind::Output) continue;
            if (connection_exists(in_node, out_node) || connection_exists(out_node, in_node)) continue;

            possible_pairs.emplace_back(in_node, out_node);
        }
    }

    if (possible_pairs.empty()) {
        return;
    }

    for (const auto& [in_node, out_node] : possible_pairs) {
        if (random_gen::generate() < adding_rate && !would_create_cycle(in_node, out_node)) {
            connections.push_back(ConnectionGene{
                in_node,
                out_node,
                1.0,
                true
            });
            dirty_topological_order = true;
        }
    }
}

void Genome::remove_connection_mutation(float removal_rate) {
    for (int i = int(connections.size()) - 1; i >= 0; i--) {
        if (random_gen::generate() < removal_rate) {
            connections.erase(connections.begin() + i);
            dirty_topological_order = true;
        }
    }
}

void Genome::add_node_mutation() {
    std::vector<int> enabled_indices;
    for (std::size_t i = 0; i < connections.size(); ++i) {
        if (connections[i].enabled) {
            enabled_indices.push_back(static_cast<int>(i));
        }
    }

    if (enabled_indices.empty()) {
        return;
    }

    int chosen_index = random_gen::random_choice(enabled_indices);
    ConnectionGene& conn = connections[chosen_index];
    conn.enabled = false;
    auto out_node = conn.out_node;
    auto weight = conn.weight;

    int new_node_id = node_id_gen.new_id();
    nodes[new_node_id] = NodeGene(
        new_node_id,
        NodeKind::Hidden,
        random_gen::random_choice(available_aggregation_funcs()),
        random_gen::random_choice(available_activation_funcs()),
        0.0
        );

    connections.push_back(ConnectionGene{
        conn.in_node,
        new_node_id,
        1.0,
        true
    });

    connections.push_back(ConnectionGene{
        new_node_id,
        out_node,
            weight,
        true
    });
    connections.erase(connections.begin() + chosen_index);
    dirty_topological_order = true;
}

void Genome::remove_node_mutation() {
    std::vector<int> node_ids;
    for (const auto& [id, _] : nodes) {
        node_ids.push_back(id);
    }

    if (node_ids.empty()) return;

    int removed_id = random_gen::random_choice(node_ids);
    if (nodes.at(removed_id).kind != NodeKind::Hidden) return;

    NodeGene removed = nodes.at(removed_id);
    nodes.erase(removed_id);

    std::vector<ConnectionGene> new_connections;
    for (const auto& conn : connections) {
        if (conn.in_node != removed.id && conn.out_node != removed.id) {
            new_connections.push_back(conn);
        }
    }
    connections = std::move(new_connections);
    dirty_topological_order = true;
}

void Genome::mutate(float weight_mutation_rate, float add_connection_rate, float add_node_rate, float remove_connection_rate, float remove_node_rate, float activation_function_rate, float max_weight_perturb) {

    mutate_weights(weight_mutation_rate, max_weight_perturb, activation_function_rate);

    if (random_gen::generate() < add_node_rate) {
        add_node_mutation();
    }
    if (random_gen::generate() < remove_node_rate) {
        remove_node_mutation();
    }
    add_connection_mutation(add_connection_rate);
    remove_connection_mutation(remove_connection_rate);

}

std::string Genome::get_symbolic_formula(const std::vector<std::string> &input_names, std::optional<std::vector<int> > output_ids) const {
    std::vector<int> out_ids = output_ids.has_value() ? *output_ids : output_node_ids;
    if (out_ids.empty()) return "";

    const NodeGene& endNode = nodes.at(out_ids[0]);
    std::string formula;

    if (std::find(input_node_ids.begin(), input_node_ids.end(), endNode.id) != input_node_ids.end()) {
        auto it = std::find(input_node_ids.begin(), input_node_ids.end(), endNode.id);
        std::size_t index = static_cast<std::size_t>(std::distance(input_node_ids.begin(), it));
        formula = input_names.at(index);
    } else {
        std::vector<ConnectionGene> conns;
        for (const auto& conn : connections) {
            if (conn.out_node == endNode.id && conn.enabled) {
                conns.push_back(conn);
            }
        }

        std::sort(conns.begin(), conns.end(),
                  [](const ConnectionGene& a, const ConnectionGene& b) {
                      return a.in_node < b.in_node;
                  });

        if (!conns.empty()) {
            std::ostringstream oss;
            for (std::size_t i = 0; i < conns.size(); ++i) {
                auto& conn = conns[i];
                if (i > 0) oss << (nodes.at(conn.out_node).aggregation == AggregationKind::Sum ? "+" : "*");
                oss << "("
                    << get_symbolic_formula(input_names, std::vector<int>{conn.in_node})
                    << ")";
                if (abs(conn.weight) - 1.0 > 0.1) {
                    oss << "*" << std::fixed << std::setprecision(1)
                    << conn.weight;
                }
            }
            formula = oss.str();
        } else {
            formula = "";
        }
    }

    std::ostringstream result;
    if (endNode.activation != ActivationFunction::Linear) result << get_activation_name(endNode.activation) << "(";
    result << formula;
    if (endNode.activation != ActivationFunction::Linear) result << ")";

    if (abs(endNode.bias) > 0.1) {
        result << (endNode.bias > 0 ? "+" : "-")
               << std::fixed << std::setprecision(1)
               << std::abs(endNode.bias);
    }
    return result.str();
}

std::string Genome::summarize() const {
    int enabled_connections = 0;
    for (const auto& c : connections) {
        if (c.enabled) ++enabled_connections;
    }

    int hidden_nodes = 0;
    for (const auto& [_, n] : nodes) {
        if (n.kind == NodeKind::Hidden) ++hidden_nodes;
    }

    std::vector<int> useful = get_useful_nodes();

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4)
        << "fitness=" << fitness << ", "
        << "nodes=" << nodes.size()
        << " (hidden=" << hidden_nodes
        << ", useful=" << useful.size() << "), "
        << "connections=" << connections.size()
        << " (enabled=" << enabled_connections << "), "
        << "Gen=" << generation;
    return oss.str();
}

Population &Population::init(int nb_inputs, int nb_outputs) {
    genomes.clear();
    for (int i = 0; i < population_size; ++i) {
        genomes.emplace_back(nb_inputs, nb_outputs);
    }
    current_best_genome = genomes[0];
    return *this;
}

std::vector<Genome> Population::reproduce() {
    std::sort(genomes.begin(), genomes.end(),
              [](const Genome& a, const Genome& b) {
                  return a.fitness > b.fitness;
              });

    int elite_count = std::max(1, static_cast<int>(genomes.size() * elite_fraction));
    std::vector<Genome> elites(genomes.begin(), genomes.begin() + elite_count);

    std::vector<Genome> new_population = elites;

    while (new_population.size() < genomes.size()) {
        const Genome& parent = random_gen::random_choice(elites);
        Genome child = parent.copy();
        child.generation = current_generation + 1;
        child.mutate(
            weight_mutation_rate,
            add_connection_rate,
            add_node_rate,
            remove_connection_rate,
            remove_node_rate,
            activation_function_rate,
            max_weight_perturb
            );
        new_population.push_back(child);
    }

    genomes = new_population;
    return new_population;
}

float sigmoid(float x) {
    return x > 10.0 ? 1.0 : (x < -10.0 ? 0.0 : 1.0 / (1.0 + std::exp(-x)));
}

float reLU(float x) {
    return std::max(0.f, x);
}

NodeGene::NodeGene(int id, NodeKind kind, AggregationKind aggregation, ActivationFunction activation, float bias)
    : id(id), kind(kind), aggregation(aggregation), activation(activation), bias(bias)
{}
