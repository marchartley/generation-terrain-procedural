from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Callable, Any, Optional, Generator
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


# ----------------------------
# Activation
# ----------------------------
def sigmoid(x: float) -> float:
    return 1 if x > 10 else 0 if x < -10 else 1.0 / (1.0 + math.exp(-x))

def reLU(x: float) -> float:
    return max(0.0, x)

def tanh(x: float) -> float:
    return math.tanh(x)

def linear(x: float) -> float:
    return x

def sin(x: float) -> float:
    return math.sin(x)

def cos(x: float) -> float:
    return math.cos(x)

# ----------------------------
# Global node ID manager
# ----------------------------

class NodeIDGenerator:
    def __init__(self, start: int = 0) -> None:
        self.next_id = start

    def new_id(self) -> int:
        value = self.next_id
        self.next_id += 1
        return value

    def copy(self):
        return NodeIDGenerator(self.next_id)

# ----------------------------
# Gene definitions
# ----------------------------

@dataclass
class NodeGene:
    id: int
    kind: str  # 'input', 'hidden', 'output'
    activation: str
    bias: float = 0.0


@dataclass
class ConnectionGene:
    in_node: int
    out_node: int
    weight: float
    enabled: bool


# ----------------------------
# Genome
# ----------------------------

# @dataclass
class Genome:

    ACTIVATION_FUNCTIONS = {
        'sigmoid': sigmoid,
        'reLU': reLU,
        'tanh': tanh,
        'linear': linear,
        # 'sin': sin,
        'cos': cos
    }
    def __init__(self, nb_inputs: int = 0, nb_outputs: int = 0) -> None:
        self.nodes: Dict[int, NodeGene] = dict()
        self.connections: List[ConnectionGene] = list()
        self.fitness: float = 0.0
        self.generation: int = 0
        self.input_node_ids: List[int] = list()
        self.output_node_ids: List[int] = list()
        self.node_id_gen: NodeIDGenerator = NodeIDGenerator(nb_inputs + nb_outputs)

        for node_id in range(nb_inputs):
            self.nodes[node_id] = NodeGene(id=node_id, kind="input", activation='linear')

        for node_id in range(nb_inputs, nb_inputs + nb_outputs):
            self.nodes[node_id] = NodeGene(id=node_id, kind="output", activation='linear')

        self.get_input_output_node_ids()

        for in_node in self.input_node_ids:
            for out_node in self.output_node_ids:
                self.connections.append(ConnectionGene(
                    in_node=in_node,
                    out_node=out_node,
                    weight=random.uniform(-1.0, 1.0),
                    enabled=True
                ))

    def get_input_output_node_ids(self) -> None:
        self.input_node_ids = [i for i, n in self.nodes.items() if n.kind == 'input']
        self.output_node_ids = [i for i, n in self.nodes.items() if n.kind == 'output']

    def copy(self) -> "Genome":
        new_genome = Genome()
        new_genome.nodes = {
            node_id: NodeGene(id=node.id, kind=node.kind, activation=node.activation, bias=node.bias)
            for node_id, node in self.nodes.items()
        }
        new_genome.input_node_ids = self.input_node_ids
        new_genome.output_node_ids = self.output_node_ids
        new_genome.node_id_gen = self.node_id_gen.copy()
        new_genome.connections = [ConnectionGene(
                in_node=conn.in_node,
                out_node=conn.out_node,
                weight=conn.weight,
                enabled=conn.enabled
            )
            for conn in self.connections
        ]
        new_genome.fitness = self.fitness
        new_genome.generation = self.generation
        return new_genome

    # ----------------------------
    # Network execution
    # ----------------------------
    def topological_order(self) -> List[int]:
        """
        Compute a topological order from enabled forward connections.

        This demo avoids recurrent connections and only allows edges
        from lower node IDs to higher node IDs, which makes this easier.
        """
        incoming_count = {node_id: 0 for node_id in self.nodes}
        outgoing: Dict[int, List[int]] = {node_id: [] for node_id in self.nodes}

        for conn in self.connections:
            if not conn.enabled:
                continue
            outgoing[conn.in_node].append(conn.out_node)
            incoming_count[conn.out_node] += 1

        queue = [node_id for node_id, count in incoming_count.items() if count == 0]
        order: List[int] = []

        while queue:
            node = queue.pop(0)
            order.append(node)
            for neighbor in outgoing[node]:
                incoming_count[neighbor] -= 1
                if incoming_count[neighbor] == 0:
                    queue.append(neighbor)

        if len(order) != len(self.nodes):
            raise ValueError("Cycle detected. This demo only supports acyclic networks.")

        return order

    def get_useful_nodes(self) -> List[Genome]:
        used_nodes = set()

        queue = set(self.output_node_ids)
        while len(queue) > 0:
            node = queue.pop()
            used_nodes.add(node)
            for conn in self.connections:
                if not conn.enabled: continue
                if conn.out_node == node:
                    queue.add(conn.in_node)
        return list(used_nodes)


    def activate_network(
        self,
        input_values: List[float]
    ) -> List[float]:
        """
        Feed-forward evaluation of a genome.
        """
        values: Dict[int, float] = {node_id: 0.0 for node_id in self.nodes}

        if len(input_values) != len(self.input_node_ids):
            raise ValueError("Input length does not match number of input nodes.")

        for node_id, value in zip(self.input_node_ids, input_values):
            values[node_id] = value

        order = self.topological_order()

        # Group incoming enabled connections by destination
        incoming_map: Dict[int, List[ConnectionGene]] = {node_id: [] for node_id in self.nodes}
        for conn in self.connections:
            if conn.enabled:
                incoming_map[conn.out_node].append(conn)

        for node_id in order:
            node = self.nodes[node_id]
            total = 0.0
            if node.kind in ("input", ):
                total = values[node_id]
            for conn in incoming_map[node_id]:
                total += values[conn.in_node] * conn.weight

            values[node_id] = self.ACTIVATION_FUNCTIONS[self.nodes[node_id].activation](total) + self.nodes[node_id].bias

        return [values[node_id] for node_id in self.output_node_ids]


    # ----------------------------
    # Mutation
    # ----------------------------

    def mutate_weights(self,
                       weight_mutation_rate: float,
                       max_weight_perturb: float,
                       activation_function_rate: float) -> None:
        for conn in self.connections:
            if random.random() < weight_mutation_rate:
                conn.weight += random.uniform(-max_weight_perturb, max_weight_perturb)

        for node in self.nodes.values():
            if random.random() < weight_mutation_rate:
                node.bias += random.uniform(-max_weight_perturb, max_weight_perturb)

            if node.kind in ("input", "output") : continue
            if random.random() < activation_function_rate:
                node.activation = random.choice(list(self.ACTIVATION_FUNCTIONS.keys()))


    def connection_exists(self,
                          in_node: int,
                          out_node: int) -> bool:
        for conn in self.connections:
            if conn.in_node == in_node and conn.out_node == out_node:
                return True
        return False

    def would_create_cycle(self,
                           in_node: int,
                           out_node: int) -> bool:
        """
        Return True if adding edge in_node -> out_node would create a cycle.

        That happens exactly when there is already a path:
            out_node -> ... -> in_node
        """
        if in_node == out_node:
            return True

        # Build adjacency list from currently enabled connections
        adjacency: Dict[int, List[int]] = {node_id: [] for node_id in self.nodes}
        for conn in self.connections:
            if conn.enabled:
                adjacency[conn.in_node].append(conn.out_node)

        # DFS/BFS from out_node looking for in_node
        stack = [out_node]
        visited = set()

        while stack:
            node = stack.pop()
            if node == in_node:
                return True
            if node in visited:
                continue
            visited.add(node)
            stack.extend(adjacency[node])

        return False

    def add_connection_mutation(self, adding_rate: float) -> None:
        node_ids = sorted(self.nodes.keys())
        possible_pairs = []

        for in_node in node_ids:
            for out_node in node_ids:
                if in_node >= out_node:
                    continue

                in_kind = self.nodes[in_node].kind
                out_kind = self.nodes[out_node].kind

                # Prevent connections INTO input/bias nodes
                if out_kind in ("input", ):
                    continue

                # Prevent connections FROM output nodes in this toy demo
                if in_kind == "output":
                    continue

                if self.connection_exists(in_node, out_node):
                    continue

                possible_pairs.append((in_node, out_node))

        if not possible_pairs:
            return

        for new_pair in possible_pairs:
            in_node, out_node = new_pair
            if random.random() < adding_rate and not self.would_create_cycle(in_node, out_node):
                self.connections.append(ConnectionGene(
                    in_node=in_node,
                    out_node=out_node,
                    weight=random.uniform(-1.0, 1.0),
                    enabled=True
                ))

    def remove_connection_mutation(self, removal_rate: float) -> None:
        for conn in self.connections:
            # if len(self.connections) == 0:
            #     return
            # removed_index = random.choice(self.connections)
            # self.connections.pop(self.connections.index(removed_index))
            if random.random() < removal_rate:
                conn.enabled = False


    def add_node_mutation(self) -> None:
        enabled_connections = [c for c in self.connections if c.enabled]
        if not enabled_connections:
            return

        conn = random.choice(enabled_connections)
        conn.enabled = False

        new_node_id = self.node_id_gen.new_id()
        self.nodes[new_node_id] = NodeGene(id=new_node_id, kind="hidden", activation=random.choice(list(self.ACTIVATION_FUNCTIONS.keys())))

        self.connections.append(ConnectionGene(
            in_node=conn.in_node,
            out_node=new_node_id,
            weight=1.0,
            enabled=True
        ))
        self.connections.append(ConnectionGene(
            in_node=new_node_id,
            out_node=conn.out_node,
            weight=conn.weight,
            enabled=True,
        ))


    def remove_node_mutation(self) -> None:
        removed_id = random.choice(list(self.nodes.keys()))
        if self.nodes[removed_id].kind != 'hidden': return
        removed = self.nodes.pop(removed_id)
        new_connections = self.connections.copy()
        for i, conn in reversed(list(enumerate(self.connections))):
            if conn.in_node == removed.id or conn.out_node == removed.id:
                new_connections.pop(i)
        self.connections = new_connections


    def mutate(self,
                weight_mutation_rate: float,
                add_connection_rate: float,
                add_node_rate: float,
                remove_connection_rate: float,
                remove_node_rate: float,
                activation_function_rate: float,
                max_weight_perturb: float
    ) -> None:
        self.mutate_weights(weight_mutation_rate, max_weight_perturb, activation_function_rate)

        if random.random() < 1.0: #add_connection_rate:
            self.add_connection_mutation(add_connection_rate)

        if random.random() < add_node_rate:
            self.add_node_mutation()

        if random.random() < 1.0: # remove_connection_rate:
            self.remove_connection_mutation(remove_connection_rate)

        if random.random() < remove_node_rate:
            self.remove_node_mutation()

    def get_symbolic_formula(self,
            input_names: List[str],
            output_node_ids: List[int] = None) -> str:

        formula = ""
        if output_node_ids is None: output_node_ids = self.output_node_ids
        endNode = self.nodes[output_node_ids[0]]

        if endNode.id in self.input_node_ids:
            formula = input_names[self.input_node_ids.index(endNode.id)]
        else:
            connections: List[ConnectionGene] = []
            for conn in self.connections:
                if conn.out_node == endNode.id and conn.enabled:
                    connections.append(conn)
            connections = sorted(connections, key=lambda c: c.in_node)
            if len(connections) > 0:
                formula = "+".join(
                    [f"({self.get_symbolic_formula(input_names, [conn.in_node])})*{conn.weight:.1f}" for conn in
                     connections])
            else:
                formula = ""
        formula = f"{endNode.activation + '(' if endNode.activation != 'linear' else ''}{formula}{')' if endNode.activation != 'linear' else ''}{'+' if endNode.bias > 0 else '-'}{abs(endNode.bias):.1f}"
        return formula

    def summarize(self) -> str:
        enabled_connections = sum(1 for c in self.connections if c.enabled)
        hidden_nodes = sum(1 for n in self.nodes.values() if n.kind == "hidden")

        useful = self.get_useful_nodes()
        return (
            f"fitness={self.fitness:.4f}, "
            f"nodes={len(self.nodes)} (hidden={hidden_nodes}, useful={len(useful)}), "
            f"connections={len(self.connections)} "
            f"(enabled={enabled_connections}), "
            f"Gen={self.generation} "
        )

@dataclass
class Population:
    population_size : int = 50
    generations : int = 10

    elite_fraction : float = 0.1
    weight_mutation_rate : float = 0.8
    add_connection_rate : float = 0.5
    add_node_rate : float = 0.5
    remove_connection_rate : float = 0.1
    remove_node_rate : float = 0.3
    activation_function_rate : float = 0.2
    network_size_loss_factor : float = 1.01

    max_weight_perturb : float = 0.5

    current_generation: int = 0

    genomes : List[Genome] = field(default_factory=list)
    current_best_genome : Genome = None

    def init(self, nb_inputs: int, nb_outputs: int) -> 'Population':
        self.genomes = [Genome(nb_inputs, nb_outputs) for _ in range(self.population_size)]
        self.current_best_genome = self.genomes[0]
        return self

    def reproduce(self) -> List[Genome]:
        population = sorted(self.genomes, key=lambda g: g.fitness, reverse=True)

        elite_count = max(1, int(len(population) * self.elite_fraction))
        elites = population[:elite_count]

        new_population: List[Genome] = [g for g in elites]

        while len(new_population) < len(population):
            parent = random.choice(elites)
            child = parent.copy()
            child.generation = self.current_generation + 1
            child.mutate(self.weight_mutation_rate, self.add_connection_rate, self.add_node_rate, self.remove_connection_rate,
                         self.remove_node_rate, self.activation_function_rate, self.max_weight_perturb)
            new_population.append(child)

        self.genomes = new_population
        return new_population

    def run(self, fitness_function: Callable[[Genome], float]) -> Generator[Genome, Any, None]:
        self.current_best_genome = self.genomes[0]
        for generation in range(self.generations):
            self.current_generation += 1
            # rebuild_dataset(generation)
            for genome in self.genomes:
                genome.fitness = fitness_function(genome)
                # evaluate_fitness(genome, TARGET_POSITIONS, NETWORK_SIZE_LOSS_FACTOR)

            population = sorted(self.genomes, key=lambda g: g.fitness, reverse=True)
            best = population[0]

            if self.current_best_genome is None or best.fitness > self.current_best_genome.fitness:
                self.current_best_genome = best.copy()
            yield self.current_best_genome
            self.reproduce()


def create_inputs_game(
    pos: np.ndarray,
    vel: np.ndarray,
    current_target: int,
    target_positions: List[np.ndarray]
) -> List[float] :
    if current_target + 1 < len(target_positions):
        next_target = target_positions[current_target + 1]
    else:
        next_target = target_positions[current_target] * 2.0 - target_positions[current_target - 1]
    return [pos[0], pos[1], vel[0], vel[1], target_positions[current_target][0], target_positions[current_target][1], next_target[0], next_target[1]]

def play_game(genome: Genome, target_positions: List[np.ndarray], max_time = 100) -> tuple[list[Any], list[Any], int, int, float]:
    p = np.array(target_positions[0])
    v = np.array([0.0, 0.0])
    wind = np.array([0.0, 0.0])

    current_target = 1
    win = False
    time_left = max_time

    Xs = []
    Ys = []
    dt = 0.5
    max_accel = 5.0
    max_speed = 3.0
    score = 0.0

    while not win and time_left > 0:
        Xs += [p[0]]
        Ys += [p[1]]
        inp = create_inputs_game(p, v, current_target, target_positions)
        predictions = genome.activate_network(inp)
        accel = np.array([predictions[0], predictions[1]])
        if np.linalg.norm(accel) > max_accel:
            accel = accel * max_accel / float(np.linalg.norm(accel))
        v += (accel + wind) * dt
        if np.linalg.norm(v) > max_speed:
            v = (v * max_speed) / float(np.linalg.norm(v))
        # v = accel
        p += v * dt
        time_left -= 1

        if np.linalg.norm(p - target_positions[current_target]) < 2.0:
            current_target += 1
            score += 100
            if current_target >= len(target_positions): win = True
    score += time_left
    if current_target < len(target_positions):
        score -= np.linalg.norm(target_positions[current_target] - p)

    return Xs, Ys, current_target, time_left, score

# ----------------------------
# Fitness evaluation:
# ----------------------------

def evaluate_fitness(
    genome: Genome,
    target_positions: List[np.ndarray],
    network_size_loss:float = 1.0
) -> float:
    _, _, current_target, time_left, score = play_game(genome, target_positions)
    fitness = score # time_left + current_target * 100
    genome.fitness = fitness
    return genome.fitness

    MSE = 0.0
    for inputs, expected in DATASET:
        predicted = genome.activate_network(inputs)
        MSE += sum(abs(p - t)**2 for p, t in zip(predicted, expected))

    MSE /= len(DATASET)
    # Higher is better
    genome.fitness = -MSE * NETWORK_SIZE_LOSS_FACTOR**(len(genome.nodes)/(len(genome.input_node_ids + genome.output_node_ids)))
    return genome.fitness

# ----------------------------
# Evolution loop
# ----------------------------


def summarize_genome(genome: Genome) -> str:
    enabled_connections = sum(1 for c in genome.connections if c.enabled)
    hidden_nodes = sum(1 for n in genome.nodes.values() if n.kind == "hidden")

    useful = genome.get_useful_nodes()
    return (
        f"fitness={genome.fitness:.4f}, "
        f"nodes={len(genome.nodes)} (hidden={hidden_nodes}, useful={len(useful)}), "
        f"connections={len(genome.connections)} "
        f"(enabled={enabled_connections}), "
        f"Gen={genome.generation} "
    )


def draw_all_predictions (
    population: List[Genome],
    target_positions: List[np.ndarray],
    error_history: List[float],
    axis: plt.Axes,
    error_axis: plt.Axes
) -> None:
    for genome in population:
        Xs, Ys, current_target, time_left, score = play_game(genome, target_positions)
        axis.plot(Xs, Ys, label = "Player", alpha = 0.1, color ="gray")

def print_predictions(
    genome: Genome,
    target_positions: List[np.ndarray],
    error_history: List[float],
    axis: plt.Axes,
    error_axis: plt.Axes
) -> None:
    Xs, Ys, current_target, time_left, score = play_game(genome, target_positions, max_time=100)
    axis.plot(Xs, Ys, label = "Player", alpha = 0.5, color ="red")
    target_colors = ["green" if i < current_target else "red" for i in range(len(target_positions))]
    axis.scatter([t[0] for t in target_positions], [t[1] for t in target_positions], label="Targets", color=target_colors)
    for i, t in enumerate(target_positions, 1):
        axis.text(t[0], t[1], str(i))
    axis.set_xlim(min([t[0] for t in target_positions] + Xs) - 2.0, max([t[0] for t in target_positions] + Xs) + 2.0)
    axis.set_ylim(min([t[1] for t in target_positions] + Ys) - 2.0, max([t[1] for t in target_positions] + Ys) + 2.0)


    error_history.append(max(score, 0.001))

    error_axis.plot(error_history[-200:])
    # error_axis.set_yscale("log")
    return




    print("\nBest network predictions on dataset:")
    MSE = 0
    GT_x = []
    GT_y = []
    Ys = []
    INs = []
    for inputs, expected in DATASET:
        predicted = activate_network(
            genome,
            inputs,
            input_node_ids,
            output_node_ids,
        )[0]
        GT_x.append(inputs[0])
        GT_y.append(expected[0])
        # Ys.append(predicted)
        MSE += np.abs(predicted - expected[0]).sum()

    for x in np.linspace(min(GT_x), max(GT_x), 100):
        _x = x
        predicted = activate_network(
            genome,
            create_inputs(_x),
            input_node_ids,
            output_node_ids,
        )[0]
        INs.append(_x)
        Ys.append(predicted)
    MSE /= len(DATASET)
    error_history.append(MSE)

    error_axis.plot(error_history)
    error_axis.set_yscale("log")
    axis.plot(INs, Ys, label="Predicted", color="orange")
    print(f"MSE: {MSE:.4f}")
    print(f"Approx. : {get_symbolic_formula(genome, input_node_ids, output_node_ids)}")
    axis.scatter(GT_x, GT_y, label="Ground truth")
    axis.set_xlabel("Expectation")
    axis.set_ylabel("Prediction")
    axis.legend()
    # test_inputs = create_inputs(6.5)
    # check_value = activate_network(genome, test_inputs, input_node_ids, output_node_ids)
    # print(f"check_value={check_value}")
    # double_check = activate_network(genome, test_inputs, input_node_ids, output_node_ids)

def visualize_genome(genome: Genome, axis: plt.Axes, input_names: List[str]) -> None:
    G = nx.DiGraph()

    node_labels = {}
    # Add nodes with type
    for node_id, node in genome.nodes.items():
        G.add_node(node_id, kind=node.kind)
        node_labels[node_id] = (input_names[node_id] if node.kind == "input" else node.activation) + f"\n{'+' if node.bias >= 0 else ''}{node.bias:.2f}"

    # Add edges (only enabled ones)
    edge_cols = []
    edge_widths = []
    for conn in genome.connections:
        if conn.enabled:
            G.add_edge(
                conn.in_node,
                conn.out_node,
                weight=round(conn.weight, 2)
            )
            edge_cols.append("red" if conn.weight < 0 else "green")
            edge_widths.append(abs(conn.weight))

    # ----------------------------
    # Layout (important!)
    # ----------------------------
    # Separate nodes by type
    input_nodes = [n for n in genome.nodes if genome.nodes[n].kind == "input"]
    hidden_nodes = [n for n in genome.nodes if genome.nodes[n].kind == "hidden"]
    output_nodes = [n for n in genome.nodes if genome.nodes[n].kind == "output"]

    # Fixed positions only for input/output
    fixed_pos = {}

    for i, node in enumerate(input_nodes):
        fixed_pos[node] = (-1, -i)

    for i, node in enumerate(output_nodes):
        fixed_pos[node] = (1, -i)

    # Compute layout:
    # - give known positions for input/output
    # - fix them
    # - hidden nodes are placed automatically
    pos = nx.spring_layout(
        G,
        pos=fixed_pos,
        fixed=input_nodes + output_nodes
    )

    useful_nodes = genome.get_useful_nodes()

    minX, maxX = [1000, -1000]
    minY, maxY = [1000, -1000]
    for n, p in pos.items():
        if n not in (input_nodes + output_nodes):
            minX, maxX = min(minX, p[0]), max(maxX, p[0])
            minY, maxY = min(minY, p[1]), max(maxY, p[1])

    if len(hidden_nodes) <= 1 :
        for n, p in pos.items():
            if n not in (input_nodes + output_nodes):
                pos[n][0] = 0.0
                pos[n][1] = 0.0
    else:
        for n, p in pos.items():
            if n not in (input_nodes + output_nodes):
                pos[n][0] = -.5 + ((p[0] - minX) / (maxX - minX)) * 1.0
                pos[n][1] = -.5 + ((p[1] - minY) / (maxY - minY)) * 1.0

    # ----------------------------
    # Node colors
    # ----------------------------
    color_map = []
    for iNode, node in enumerate(G.nodes()):
        if node not in useful_nodes:
            color_map.append("white")
            continue
        kind = genome.nodes[node].kind
        if kind == "input":
            color_map.append("lightblue")
        elif kind == "hidden":
            color_map.append("gray")
        elif kind == "output":
            color_map.append("lightgreen")

    # ----------------------------
    # Draw graph
    # ----------------------------
    # plt.figure(figsize=(8, 6))

    nx.draw(
        G,
        pos,
        with_labels=False,
        node_color=color_map,
        node_size=800,
        # arrows=True
        ax = axis,
        edge_color=edge_cols,
        width=edge_widths,
    )

    # Draw edge labels (weights)
    edge_labels = nx.get_edge_attributes(G, "weight")
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, ax=axis)
    nx.draw_networkx_labels(G, pos, node_labels, font_size=8, ax=axis)

    axis.set_title("NEAT Genome Network")
    axis.axis("off")




class NeatViewer:
    def __init__(self, input_names: List[str], population: Population) -> None:
        fig = plt.figure()
        self.brain_ax: plt.Axes = fig.add_subplot(1, 3, 1)
        self.game_ax: plt.Axes = fig.add_subplot(1, 3, 2)
        self.fitness_ax: plt.Axes = fig.add_subplot(1, 3, 3)

        self.fitness_history: List[float] = []
        self.input_names: List[str] = input_names
        self.population = population

    def show_brain(self, genome: Genome):
        self.brain_ax.clear()

        useful_nodes = genome.get_useful_nodes()


        G = nx.DiGraph()

        node_labels = {}
        # Add nodes with type
        for node_id in useful_nodes:
            node = genome.nodes[node_id]
            G.add_node(node_id, kind=node.kind)
            node_labels[node_id] = (self.input_names[node_id] if node.kind == "input" else node.activation) + f"\n{'+' if node.bias >= 0 else ''}{node.bias:.2f}"

        # Add edges (only enabled ones)
        edge_cols = []
        edge_widths = []
        for conn in genome.connections:
            if conn.enabled and conn.in_node in useful_nodes and conn.out_node in useful_nodes:
                G.add_edge(
                    conn.in_node,
                    conn.out_node,
                    weight=round(conn.weight, 2)
                )
                edge_cols.append("red" if conn.weight < 0 else "green")
                edge_widths.append(abs(conn.weight))

        # ----------------------------
        # Layout (important!)
        # ----------------------------
        # Separate nodes by type
        input_nodes = [n for n in genome.nodes if genome.nodes[n].kind == "input" and n in useful_nodes]
        hidden_nodes = [n for n in genome.nodes if genome.nodes[n].kind == "hidden" and n in useful_nodes]
        output_nodes = [n for n in genome.nodes if genome.nodes[n].kind == "output" and n in useful_nodes]

        # Fixed positions only for input/output
        fixed_pos = {}

        for i, node in enumerate(input_nodes):
            fixed_pos[node] = (-1, -i)

        for i, node in enumerate(output_nodes):
            fixed_pos[node] = (1, -i)

        # Compute layout:
        # - give known positions for input/output
        # - fix them
        # - hidden nodes are placed automatically
        pos = nx.spring_layout(
            G,
            pos=fixed_pos,
            fixed=input_nodes + output_nodes
        )

        minX, maxX = [1000, -1000]
        minY, maxY = [1000, -1000]
        for n, p in pos.items():
            if n not in (input_nodes + output_nodes):
                minX, maxX = min(minX, p[0]), max(maxX, p[0])
                minY, maxY = min(minY, p[1]), max(maxY, p[1])

        if len(hidden_nodes) <= 1:
            for n, p in pos.items():
                if n not in (input_nodes + output_nodes):
                    pos[n][0] = 0.0
                    pos[n][1] = 0.0
        else:
            for n, p in pos.items():
                if n not in (input_nodes + output_nodes):
                    pos[n][0] = -.5 + ((p[0] - minX) / (maxX - minX)) * 1.0
                    pos[n][1] = -.5 + ((p[1] - minY) / (maxY - minY)) * 1.0

        # ----------------------------
        # Node colors
        # ----------------------------
        color_map = []
        for iNode, node in enumerate(G.nodes()):
            # if node not in useful_nodes:
                # color_map.append("white")
                # continue
            kind = genome.nodes[node].kind
            if kind == "input":
                color_map.append("lightblue")
            elif kind == "hidden":
                color_map.append("gray")
            elif kind == "output":
                color_map.append("lightgreen")

        # ----------------------------
        # Draw graph
        # ----------------------------
        # plt.figure(figsize=(8, 6))


        # Draw edge labels (weights)
        edge_labels = nx.get_edge_attributes(G, "weight")
        nx.draw_networkx_labels(G, pos, node_labels, font_size=8, ax=self.brain_ax)
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, ax=self.brain_ax)
        nx.draw(
            G,
            pos,
            with_labels=False,
            node_color=color_map,
            node_size=800,
            # arrows=True
            ax=self.brain_ax,
            edge_color=edge_cols,
            width=edge_widths,
        )

        self.brain_ax.set_title("Best Network")
        self.brain_ax.axis("off")

    def show_game(self):
        pass

    def show_error(self):
        self.fitness_ax.clear()
        self.fitness_ax.plot(self.fitness_history)

    def game_loop(self):
        pass

    def __prepare_start(self):
        plt.ion()
        plt.show()
        plt.pause(0.01)

    def __prepare_end(self):
        plt.ioff()
        plt.show()

    def start(self):
        self.__prepare_start()

        self.game_loop()

        self.__prepare_end()


def main() -> None:
    # ----------------------------
    # Configuration
    # ----------------------------
    population = Population(
        population_size = 50,
        generations = 10,
        elite_fraction = 0.1,

        weight_mutation_rate = 0.8,
        add_connection_rate = 0.5,
        add_node_rate = 0.5,
        remove_connection_rate = 0.1,
        remove_node_rate = 0.3,
        activation_function_rate = 0.2,

        max_weight_perturb = 0.5
    ).init(nb_inputs=8, nb_outputs=2)

    NETWORK_SIZE_LOSS_FACTOR = 1.01


    def rebuild_dataset(generation: int = 0) -> List[Any]:
        def create_inputs(x: float) -> List[float]:
            return [x]  # , x**2, math.sqrt(abs(x))]

        def create_outputs(x: float) -> List[float]:
            return [math.cos(x) * 8.0 + 4.0]
        x0 = 0.0
        x1 = 10.0
        off = 0.0  # math.sin(generation * 0.1) * 10
        random_values = np.random.uniform(x0 + off, x1 + off, size=200).tolist()
        dataset = [(create_inputs(x), create_outputs(x + (random.random() - .5) * 1.0)) for x in random_values]
        return dataset

    # INPUT_NAMES = ['x'] #, 'x^2', '$\sqrt{|x|}$']
    INPUT_NAMES = ['$p_x$', '$p_y$', '$v_x$', '$v_y$', '$t_x$', '$t_y$', '$t_{x+1}$', '$t_{y+1}$']

    NB_TARGETS = 10
    TARGET_POSITIONS: List[np.ndarray] = [
        np.random.uniform(-20.0, 20.0, size=2)
        for _ in range(NB_TARGETS)
    ]

    error_history: List[float] = []
    # dataset = rebuild_dataset()
    plt.ion()
    fig = plt.figure(figsize=(8, 6))
    current_best_visual = fig.add_subplot(1, 3, 1)
    current_best_results = fig.add_subplot(1, 3, 2)
    current_error_ax = fig.add_subplot(1, 3, 3)
    plt.show()
    plt.pause(0.01)

    best_ever = population.genomes[0]

    for current_best in population.run(lambda genome : evaluate_fitness(genome, TARGET_POSITIONS, NETWORK_SIZE_LOSS_FACTOR)):
        best_ever = current_best
        generation = population.current_generation

        rebuild_dataset(generation)
        current_best_visual.clear()
        current_best_results.clear()
        current_error_ax.clear()
        visualize_genome(best_ever, current_best_visual, INPUT_NAMES)
        draw_all_predictions(population.genomes, TARGET_POSITIONS, error_history, current_best_results, current_error_ax)
        print_predictions(best_ever, TARGET_POSITIONS, error_history, current_best_results, current_error_ax)
        plt.pause(0.01)

        print(f"Generation {generation:02d}: {summarize_genome(best_ever)}")

    print("\nBest genome found:")
    print(summarize_genome(best_ever))

    current_best_visual.clear()
    current_best_results.clear()
    current_error_ax.clear()
    visualize_genome(best_ever, current_best_visual, INPUT_NAMES)
    best_ever.get_symbolic_formula(INPUT_NAMES)
    draw_all_predictions(population.genomes, TARGET_POSITIONS, error_history, current_best_results, current_error_ax)
    print_predictions(best_ever, TARGET_POSITIONS, error_history, current_best_results, current_error_ax)
    plt.ioff()
    plt.show()


if __name__ == "__main__":
    main()