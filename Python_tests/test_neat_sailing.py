import math
import random

from test_neat import Population, Genome, NeatViewer
import matplotlib.pyplot as plt
from typing import List, Callable, Any
import numpy as np
import networkx as nx


class NeatSailingViewer(NeatViewer):
    def __init__(self, population: Population) -> None:
        super().__init__(['$p_x$', '$p_y$', '$v_x$', '$v_y$', '$t_x$', '$t_y$', '$t_{x+1}$', '$t_{y+1}$'], population)
        self.target_positions: List[np.ndarray] = []

    def show_game(self):
        self.game_ax.clear()
        for genome in self.population.genomes:
            Xs, Ys, current_target, time_left, score = self.play_game(genome, self.target_positions)
            self.game_ax.plot(Xs, Ys, label = "Player", alpha = 0.1, color ="gray")


        Xs, Ys, current_target, time_left, score = self.play_game(self.population.current_best_genome, self.target_positions, max_time=100)
        self.game_ax.plot(Xs, Ys, label="Player", alpha=0.5, color="red")
        target_colors = ["green" if i < current_target else "red" for i in range(len(self.target_positions))]
        self.game_ax.scatter([t[0] for t in self.target_positions], [t[1] for t in self.target_positions], label="Targets",
                     color=target_colors)
        for i, t in enumerate(self.target_positions, 1):
            self.game_ax.text(t[0], t[1], str(i))
        self.game_ax.set_xlim(min([t[0] for t in self.target_positions] + Xs) - 2.0,
                      max([t[0] for t in self.target_positions] + Xs) + 2.0)
        self.game_ax.set_ylim(min([t[1] for t in self.target_positions] + Ys) - 2.0,
                      max([t[1] for t in self.target_positions] + Ys) + 2.0)

    def game_loop(self):
        NETWORK_SIZE_LOSS_FACTOR = 1.01
        NB_TARGETS = 10
        self.target_positions = [np.random.uniform(-20.0, 20.0, size=2) for _ in range(NB_TARGETS) ]
        best_ever = self.population.genomes[0]

        for current_best in self.population.run(
                lambda genome: self.evaluate_fitness(genome, self.target_positions, NETWORK_SIZE_LOSS_FACTOR)):
            best_ever = current_best
            generation = self.population.current_generation

            self.fitness_history.append(max(best_ever.fitness, 0.001))

            self.show_brain(best_ever)
            self.show_game()
            self.show_error()

            print(f"Generation {generation:02d}: {best_ever.summarize()}")
            plt.pause(0.01)

            if generation % 50 == 0:
                self.target_positions = [np.random.uniform(-20.0, 20.0, size=2) for _ in range(NB_TARGETS)]

        print("\nBest genome found:")
        print(best_ever.summarize())

    def create_inputs_game(self, pos: np.ndarray, vel: np.ndarray, current_target: int, target_positions: List[np.ndarray]) -> List[float] :
        if current_target + 1 < len(target_positions):
            next_target = target_positions[current_target + 1]
        else:
            next_target = target_positions[current_target] * 2.0 - target_positions[current_target - 1]
        return [pos[0], pos[1], vel[0], vel[1], target_positions[current_target][0], target_positions[current_target][1], next_target[0], next_target[1]]

    def play_game(self, genome: Genome, target_positions: List[np.ndarray], max_time = 100) -> tuple[list[Any], list[Any], int, int, float]:
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
            inp = self.create_inputs_game(p, v, current_target, target_positions)
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

    def evaluate_fitness(self, genome: Genome, target_positions: List[np.ndarray], network_size_loss:float = 1.0) -> float:
        _, _, current_target, time_left, score = self.play_game(genome, target_positions)
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





def main() -> None:
    population = Population(
        population_size = 100,
        generations = 1000,
        elite_fraction = 0.05,

        weight_mutation_rate = 1.0,
        add_connection_rate = 0.2,
        add_node_rate = 0.3,
        remove_connection_rate = 0.2,
        remove_node_rate = 0.3,
        activation_function_rate = 0.1,

        max_weight_perturb = 0.1
    ).init(nb_inputs=8, nb_outputs=2)

    viewer = NeatSailingViewer(population = population)
    viewer.start()


if __name__ == "__main__":
    main()