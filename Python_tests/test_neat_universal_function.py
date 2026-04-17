import math
import random
from typing import List

from matplotlib import pyplot as plt

from test_neat import NeatViewer, Population, Genome


class NeatUniversalFuncViewer(NeatViewer):
    def __init__(self, population: Population):
        super().__init__(["$x$", "$x^2$"], population)
        self.dataset: List[float] = sorted([random.random() * 10.0 for _ in range(100)])

    def make_output(self, x: float):
        return math.cos(x * 2.0) + 3

    def make_inputs(self, x: float):
        return [x, x*x]

    def show_game(self):
        self.game_ax.clear()
        self.game_ax.scatter(self.dataset, [self.make_output(x) for x in self.dataset])
        self.game_ax.plot(self.dataset, [self.population.current_best_genome.activate_network(self.make_inputs(x)) for x in self.dataset])

    def game_loop(self):
        NETWORK_SIZE_LOSS_FACTOR = 1.01
        best_ever = self.population.genomes[0]

        for current_best in self.population.run(
                lambda genome: self.evaluate_fitness(genome, NETWORK_SIZE_LOSS_FACTOR)):
            best_ever = current_best
            generation = self.population.current_generation

            self.fitness_history.append(max(best_ever.fitness, 0.001))

            self.show_brain(best_ever)
            self.show_game()
            self.show_error()

            print(f"Generation {generation:02d}: {best_ever.summarize()}")
            plt.pause(0.01)

        print("\nBest genome found:")
        print(best_ever.summarize())

    def evaluate_fitness(self, genome: Genome, network_size_loss:float = 1.0) -> float:
        error = 0.0
        for x in self.dataset:
            error += (genome.activate_network(self.make_inputs(x))[0] - self.make_output(x))**2.0
        # error /= float(len(self.dataset))
        # error /= network_size_loss**(len(genome.nodes))
        return 100.0 - error


if __name__ == "__main__":
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
    ).init(2, 1)
    # .init(nb_inputs=8, nb_outputs=2))

    # viewer = NeatSailingViewer(population = population)
    viewer = NeatUniversalFuncViewer(population)
    viewer.start()