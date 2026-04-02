from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Callable, Any
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------
# Configuration
# ----------------------------

POPULATION_SIZE = 50
GENERATIONS = 100
ELITE_FRACTION = 0.1

WEIGHT_MUTATION_RATE = 0.8
ADD_CONNECTION_RATE = 0.5
ADD_NODE_RATE = 0.5
REMOVE_CONNECTION_RATE = 0.1
REMOVE_NODE_RATE = 0.3
ACTIVATION_FUNCTION_RATE = 0.2
NETWORK_SIZE_LOSS_FACTOR = 1.01

MAX_WEIGHT_PERTURB = 0.5
RANDOM_SEED = 7

#
# DATASET = [
#     ([0.0, 0.0], [0.0]),
#     ([0.0, 1.0], [1.0]),
#     ([1.0, 0.0], [1.0]),
#     ([1.0, 1.0], [0.0]),
# ]
def create_inputs(x:float) -> List[float]:
    return [x] #, x**2, math.sqrt(abs(x))]
def create_outputs(x: float) -> List[float]:
    return [math.cos(x) * 8.0 + 4.0]

def rebuild_dataset(generation: int = 0) -> List[Any]:
    global DATASET
    x0 = 0.0
    x1 = 10.0
    off = 0.0 # math.sin(generation * 0.1) * 10
    DATASET = [(create_inputs(x), create_outputs(x + (random.random() -.5) * 1.0)) for x in np.random.uniform(x0 + off, x1 + off, size=200)]
    return DATASET

DATASET = rebuild_dataset()
# INPUT_NAMES = ['x'] #, 'x^2', '$\sqrt{|x|}$']
INPUT_NAMES = ['$p_x$', '$p_y$', '$v_x$', '$v_y$', '$t_x$', '$t_y$', '$t_{x+1}$', '$t_{y+1}$']

NB_TARGETS = 10
TARGET_POSITIONS =  [
    # np.array([0.0,  5.0]),
    # np.array([20.0, 10.0]),
    # np.array([30.0, 0.0]),
    # np.array([10.0, 30.0]),
    # np.array([30.0, -10.0])
    np.random.uniform(-20.0, 20.0, size=2)
    for _ in range(NB_TARGETS)
]



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

ACTIVATION_FUNCTIONS = {
    'sigmoid': sigmoid,
    'reLU': reLU,
    'tanh': tanh,
    'linear': linear,
    # 'sin': sin,
    'cos': cos
}


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

@dataclass
class Genome:
    nodes: Dict[int, NodeGene] = field(default_factory=dict)
    connections: List[ConnectionGene] = field(default_factory=list)
    fitness: float = 0.0
    generation: int = 0

    def copy(self) -> "Genome":
        new_genome = Genome()
        new_genome.nodes = {
            node_id: NodeGene(id=node.id, kind=node.kind, activation=node.activation, bias=node.bias)
            for node_id, node in self.nodes.items()
        }
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
# Global node ID manager
# ----------------------------

class NodeIDGenerator:
    def __init__(self, start: int = 0) -> None:
        self.next_id = start

    def new_id(self) -> int:
        value = self.next_id
        self.next_id += 1
        return value


# ----------------------------
# Network execution
# ----------------------------

def topological_order(genome: Genome) -> List[int]:
    """
    Compute a topological order from enabled forward connections.

    This demo avoids recurrent connections and only allows edges
    from lower node IDs to higher node IDs, which makes this easier.
    """
    incoming_count = {node_id: 0 for node_id in genome.nodes}
    outgoing: Dict[int, List[int]] = {node_id: [] for node_id in genome.nodes}

    for conn in genome.connections:
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

    if len(order) != len(genome.nodes):
        raise ValueError("Cycle detected. This demo only supports acyclic networks.")

    return order


def activate_network(
    genome: Genome,
    input_values: List[float],
    input_node_ids: List[int],
    output_node_ids: List[int],
) -> List[float]:
    """
    Feed-forward evaluation of a genome.
    """
    values: Dict[int, float] = {node_id: 0.0 for node_id in genome.nodes}

    if len(input_values) != len(input_node_ids):
        raise ValueError("Input length does not match number of input nodes.")

    for node_id, value in zip(input_node_ids, input_values):
        values[node_id] = value

    order = topological_order(genome)

    # Group incoming enabled connections by destination
    incoming_map: Dict[int, List[ConnectionGene]] = {node_id: [] for node_id in genome.nodes}
    for conn in genome.connections:
        if conn.enabled:
            incoming_map[conn.out_node].append(conn)

    for node_id in order:
        node = genome.nodes[node_id]
        total = 0.0
        if node.kind in ("input", ):
            total = values[node_id]
        for conn in incoming_map[node_id]:
            total += values[conn.in_node] * conn.weight

        values[node_id] = ACTIVATION_FUNCTIONS[genome.nodes[node_id].activation](total) + genome.nodes[node_id].bias

    return [values[node_id] for node_id in output_node_ids]


# ----------------------------
# Mutation
# ----------------------------

def mutate_weights(genome: Genome) -> None:
    for conn in genome.connections:
        if random.random() < WEIGHT_MUTATION_RATE:
            conn.weight += random.uniform(-MAX_WEIGHT_PERTURB, MAX_WEIGHT_PERTURB)

    for node in genome.nodes.values():
        if random.random() < WEIGHT_MUTATION_RATE:
            node.bias += random.uniform(-MAX_WEIGHT_PERTURB, MAX_WEIGHT_PERTURB)

        if node.kind in ("input", "output") : continue
        if random.random() < ACTIVATION_FUNCTION_RATE:
            node.activation = random.choice(list(ACTIVATION_FUNCTIONS.keys()))


def connection_exists(genome: Genome, in_node: int, out_node: int) -> bool:
    for conn in genome.connections:
        if conn.in_node == in_node and conn.out_node == out_node:
            return True
    return False

def would_create_cycle(genome: Genome, in_node: int, out_node: int) -> bool:
    """
    Return True if adding edge in_node -> out_node would create a cycle.

    That happens exactly when there is already a path:
        out_node -> ... -> in_node
    """
    if in_node == out_node:
        return True

    # Build adjacency list from currently enabled connections
    adjacency: Dict[int, List[int]] = {node_id: [] for node_id in genome.nodes}
    for conn in genome.connections:
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

def add_connection_mutation(genome: Genome) -> None:
    node_ids = sorted(genome.nodes.keys())
    possible_pairs = []

    for in_node in node_ids:
        for out_node in node_ids:
            if in_node >= out_node:
                continue

            in_kind = genome.nodes[in_node].kind
            out_kind = genome.nodes[out_node].kind

            # Prevent connections INTO input/bias nodes
            if out_kind in ("input", ):
                continue

            # Prevent connections FROM output nodes in this toy demo
            if in_kind == "output":
                continue

            if connection_exists(genome, in_node, out_node):
                continue

            if would_create_cycle(genome, in_node, out_node):
                continue

            possible_pairs.append((in_node, out_node))

    if not possible_pairs:
        return

    in_node, out_node = random.choice(possible_pairs)
    genome.connections.append(ConnectionGene(
        in_node=in_node,
        out_node=out_node,
        weight=random.uniform(-1.0, 1.0),
        enabled=True
    ))

def remove_connection_mutation(genome: Genome) -> None:
    if len(genome.connections) == 0:
        return
    removed_index = random.choice(genome.connections)
    genome.connections.pop(genome.connections.index(removed_index))


def add_node_mutation(
    genome: Genome,
    node_ids: NodeIDGenerator,
) -> None:
    enabled_connections = [c for c in genome.connections if c.enabled]
    if not enabled_connections:
        return

    conn = random.choice(enabled_connections)
    conn.enabled = False

    new_node_id = node_ids.new_id()
    genome.nodes[new_node_id] = NodeGene(id=new_node_id, kind="hidden", activation=random.choice(list(ACTIVATION_FUNCTIONS.keys())))

    genome.connections.append(ConnectionGene(
        in_node=conn.in_node,
        out_node=new_node_id,
        weight=1.0,
        enabled=True
    ))
    genome.connections.append(ConnectionGene(
        in_node=new_node_id,
        out_node=conn.out_node,
        weight=conn.weight,
        enabled=True,
    ))


def remove_node_mutation(
    genome: Genome
) -> None:
    removed_id = random.choice(list(genome.nodes.keys()))
    if genome.nodes[removed_id].kind != 'hidden': return
    removed = genome.nodes.pop(removed_id)
    new_connections = genome.connections.copy()
    for i, conn in reversed(list(enumerate(genome.connections))):
        if conn.in_node == removed.id or conn.out_node == removed.id:
            new_connections.pop(i)
    genome.connections = new_connections


def mutate(
    genome: Genome,
    node_ids: NodeIDGenerator,
) -> None:
    mutate_weights(genome)

    if random.random() < ADD_CONNECTION_RATE:
        add_connection_mutation(genome)

    if random.random() < ADD_NODE_RATE:
        add_node_mutation(genome, node_ids)

    if random.random() < REMOVE_CONNECTION_RATE:
        remove_connection_mutation(genome)

    if random.random() < REMOVE_NODE_RATE:
        remove_node_mutation(genome)

def create_inputs_game(
    pos: np.ndarray,
    vel: np.ndarray,
    current_target: int
) -> List[float] :
    if current_target + 1 < len(TARGET_POSITIONS):
        next_target = TARGET_POSITIONS[current_target + 1]
    else:
        next_target = TARGET_POSITIONS[current_target] * 2.0 - TARGET_POSITIONS[current_target - 1]
    return [pos[0], pos[1], vel[0], vel[1], TARGET_POSITIONS[current_target][0], TARGET_POSITIONS[current_target][1], next_target[0], next_target[1]]

def play_game(genome: Genome, max_time = 100) -> tuple[list[Any], list[Any], int, int, float]:

    p = np.array(TARGET_POSITIONS[0])
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
        inp = create_inputs_game(p, v, current_target)
        predictions = activate_network(genome, inp, [n.id for n in genome.nodes.values() if n.kind == 'input'], [n.id for n in genome.nodes.values() if n.kind == 'output'])
        accel = np.array([predictions[0], predictions[1]])
        if np.linalg.norm(accel) > max_accel:
            accel = accel * max_accel / np.linalg.norm(accel)
        v += (accel + wind) * dt
        if np.linalg.norm(v) > max_speed:
            v = (v * max_speed) / np.linalg.norm(v)
        # v = accel
        p += v * dt
        time_left -= 1

        if np.linalg.norm(p - TARGET_POSITIONS[current_target]) < 2.0:
            current_target += 1
            score += 100
            if current_target >= len(TARGET_POSITIONS): win = True
    score += time_left
    if current_target < len(TARGET_POSITIONS):
        score -= np.linalg.norm(TARGET_POSITIONS[current_target] - p)

    return Xs, Ys, current_target, time_left, score

# ----------------------------
# Fitness evaluation:
# ----------------------------

def evaluate_fitness(
    genome: Genome,
    input_node_ids: List[int],
    output_node_ids: List[int],
) -> float:
    _, _, current_target, time_left, score = play_game(genome)
    fitness = score # time_left + current_target * 100
    genome.fitness = fitness
    return genome.fitness

    MSE = 0.0
    for inputs, expected in DATASET:
        predicted = activate_network(
            genome=genome,
            input_values=inputs,
            input_node_ids=input_node_ids,
            output_node_ids=output_node_ids,
        )
        MSE += sum(abs(p - t)**2 for p, t in zip(predicted, expected))

    MSE /= len(DATASET)
    # Higher is better
    genome.fitness = -MSE * NETWORK_SIZE_LOSS_FACTOR**(len(genome.nodes)/(len(input_node_ids + output_node_ids)))
    return genome.fitness


# ----------------------------
# Population initialization
# ----------------------------

def make_initial_genome(
    input_node_ids: List[int],
    output_node_ids: List[int],
) -> Genome:
    genome = Genome()

    for node_id in input_node_ids:
        genome.nodes[node_id] = NodeGene(id=node_id, kind="input", activation='linear')

    for node_id in output_node_ids:
        genome.nodes[node_id] = NodeGene(id=node_id, kind="output", activation='linear')

    for in_node in input_node_ids:
        for out_node in output_node_ids:
            genome.connections.append(ConnectionGene(
                in_node=in_node,
                out_node=out_node,
                weight=random.uniform(-1.0, 1.0),
                enabled=True
            ))

    return genome


# ----------------------------
# Evolution loop
# ----------------------------

def reproduce(
    population: List[Genome],
    node_ids: NodeIDGenerator,
    current_generation: int,
) -> List[Genome]:
    population = sorted(population, key=lambda g: g.fitness, reverse=True)

    elite_count = max(1, int(len(population) * ELITE_FRACTION))
    elites = population[:elite_count]

    new_population: List[Genome] = [g for g in elites]

    while len(new_population) < len(population):
        parent = random.choice(elites)
        child = parent.copy()
        child.generation = current_generation + 1
        mutate(child, node_ids)
        new_population.append(child)

    return new_population

def get_symbolic_formula(
    genome: Genome,
    input_node_ids: List[int],
    output_node_ids: List[int] ) -> str:

    formula = ""
    endNode = genome.nodes[output_node_ids[0]]

    if endNode.id in input_node_ids:
        formula = INPUT_NAMES[input_node_ids.index(endNode.id)]
    else:
        connections: List[ConnectionGene] = []
        for conn in genome.connections:
            if conn.out_node == endNode.id and conn.enabled:
                connections.append(conn)
        connections = sorted(connections, key=lambda c: c.in_node)
        if len(connections) > 0:
            formula = "+".join([f"({get_symbolic_formula(genome, input_node_ids, [conn.in_node])})*{conn.weight:.1f}" for conn in connections])
        else:
            formula = ""
    formula = f"{endNode.activation + '(' if endNode.activation != 'linear' else ''}{formula}{')' if endNode.activation != 'linear' else ''}{'+' if endNode.bias > 0 else '-'}{abs(endNode.bias):.1f}"
    return formula

def get_useful_nodes(genome: Genome, output_node_ids: List[int]) -> List[Genome]:
    used_nodes = set()

    queue = set(output_node_ids)
    while len(queue) > 0:
        node = queue.pop()
        used_nodes.add(node)
        for conn in genome.connections:
            if not conn.enabled: continue
            if conn.out_node == node:
                queue.add(conn.in_node)
    return list(used_nodes)

def summarize_genome(genome: Genome) -> str:
    enabled_connections = sum(1 for c in genome.connections if c.enabled)
    hidden_nodes = sum(1 for n in genome.nodes.values() if n.kind == "hidden")

    useful = get_useful_nodes(genome, [n.id for k, n in genome.nodes.items() if n.kind == 'output'])
    return (
        f"fitness={genome.fitness:.4f}, "
        f"nodes={len(genome.nodes)} (hidden={hidden_nodes}, useful={len(useful)}), "
        f"connections={len(genome.connections)} "
        f"(enabled={enabled_connections}), "
        f"Gen={genome.generation} "
    )


def draw_all_predictions (
    population: List[Genome],
    input_node_ids: List[int],
    output_node_ids: List[int],
    error_history: List[float],
    axis: plt.Axes,
    error_axis: plt.Axis
) -> None:
    for genome in population:
        Xs, Ys, current_target, time_left, score = play_game(genome)
        axis.plot(Xs, Ys, label = "Player", alpha = 0.1, color ="gray")

def print_predictions(
    genome: Genome,
    input_node_ids: List[int],
    output_node_ids: List[int],
    error_history: List[float],
    axis: plt.Axes,
    error_axis: plt.Axis
) -> None:
    Xs, Ys, current_target, time_left, score = play_game(genome, max_time=100)
    axis.plot(Xs, Ys, label = "Player", alpha = 0.5, color ="red")
    target_colors = ["green" if i < current_target else "red" for i in range(len(TARGET_POSITIONS))]
    axis.scatter([t[0] for t in TARGET_POSITIONS], [t[1] for t in TARGET_POSITIONS], label="Targets", color=target_colors)
    for i, t in enumerate(TARGET_POSITIONS, 1):
        axis.text(t[0], t[1], str(i))
    axis.set_xlim(min([t[0] for t in TARGET_POSITIONS] + Xs) - 2.0, max([t[0] for t in TARGET_POSITIONS] + Xs) + 2.0)
    axis.set_ylim(min([t[1] for t in TARGET_POSITIONS] + Ys) - 2.0, max([t[1] for t in TARGET_POSITIONS] + Ys) + 2.0)


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

def visualize_genome(genome: Genome, axis: plt.Axes) -> None:
    G = nx.DiGraph()

    node_labels = {}
    # Add nodes with type
    for node_id, node in genome.nodes.items():
        G.add_node(node_id, kind=node.kind)
        node_labels[node_id] = (INPUT_NAMES[node_id] if node.kind == "input" else node.activation) + f"\n{'+' if node.bias >= 0 else ''}{node.bias:.2f}"

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

    useful_nodes = get_useful_nodes(genome, output_nodes)

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
    for node in G.nodes():
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

def main() -> None:
    # random.seed(RANDOM_SEED)
    error_history: List[float] = []
    rebuild_dataset()
    plt.ion()
    fig = plt.figure(figsize=(8, 6))
    current_best_visual = fig.add_subplot(1, 3, 1)
    current_best_results = fig.add_subplot(1, 3, 2)
    current_error_ax = fig.add_subplot(1, 3, 3)
    plt.show()
    plt.pause(0.01)

    # Fixed initial nodes
    node_ids = NodeIDGenerator(start=0)
    # input_node_ids = [node_ids.new_id() for _ in range(len(DATASET[0][0]))]   # x1, x2
    input_node_ids = [node_ids.new_id() for _ in range(8)]
    output_node_ids = [node_ids.new_id() for _ in range(2)]                     # y

    population = [
        make_initial_genome(input_node_ids, output_node_ids)
        for _ in range(POPULATION_SIZE)
    ]

    best_ever: Genome | None = None

    for generation in range(GENERATIONS):
        rebuild_dataset(generation)
        for genome in population:
            evaluate_fitness(genome, input_node_ids, output_node_ids)

        population.sort(key=lambda g: g.fitness, reverse=True)
        best = population[0]

        if best_ever is None or best.fitness > best_ever.fitness:
            best_ever = best.copy()
        current_best_visual.clear()
        current_best_results.clear()
        current_error_ax.clear()
        visualize_genome(best_ever, current_best_visual)
        draw_all_predictions(population, input_node_ids, output_node_ids, error_history, current_best_results, current_error_ax)
        print_predictions(best_ever, input_node_ids, output_node_ids, error_history, current_best_results, current_error_ax)
        plt.pause(0.01)

        print(f"Generation {generation:02d}: {summarize_genome(best)}")

        population = reproduce(population, node_ids, generation)

    assert best_ever is not None
    print("\nBest genome found:")
    print(summarize_genome(best_ever))

    current_best_visual.clear()
    current_best_results.clear()
    current_error_ax.clear()
    visualize_genome(best_ever, current_best_visual)
    get_symbolic_formula(best_ever, input_node_ids, output_node_ids)
    draw_all_predictions(population, input_node_ids, output_node_ids, error_history, current_best_results, current_error_ax)
    print_predictions(best_ever, input_node_ids, output_node_ids, error_history, current_best_results, current_error_ax)
    plt.ioff()
    plt.show()


if __name__ == "__main__":
    main()