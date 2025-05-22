import random
import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from concurrent.futures import ThreadPoolExecutor
import copy
import multiprocessing as mp


# Classe pour une masse
class Mass:
    def __init__(self, position, mass=10.0):
        self.position = np.array(position, dtype=float)
        self.mass = mass
        self.velocity = np.zeros(2)
        self.force = np.zeros(2)
        self.prev_position = np.array(position, dtype=float)
        self.v_half = np.zeros(2)
        self.fixed = False

# Classe pour un ressort
class Spring:
    def __init__(self, mass1_idx, mass2_idx, rest_length, stiffness):
        self.mass1_idx = mass1_idx
        self.mass2_idx = mass2_idx
        self.rest_length = rest_length
        self.stiffness = stiffness

# Initialisation de la grille
def create_mass_spring_grid(N, M, k):
    masses = []
    springs = []

    # Crée les masses avec position initiale en grille régulière
    rand_range = 0.5
    for i in range(N):
        for j in range(M):
            fixed = False
            if j == M-1: # Fix top row
                position = [i, j]
                fixed = True
            elif j > 0:
                position = [i, j]
            else:
                position = [i + (random.random() - 0.5) * rand_range, j + (random.random() - 0.5) * rand_range]
            masses.append(Mass(position))
            masses[-1].fixed = fixed

    # Crée les ressorts
    def index(i, j):
        return i * M + j

    for i in range(N):
        for j in range(M):
            current_idx = index(i, j)

            # Voisin à droite
            if j + 1 < M:
                springs.append(Spring(current_idx, index(i, j+1), 1.0, k))
            # Voisin en bas
            if i + 1 < N:
                springs.append(Spring(current_idx, index(i+1, j), 1.0, k))
            # Diagonale bas-droite
            if i + 1 < N and j + 1 < M:
                springs.append(Spring(current_idx, index(i+1, j+1), sqrt(2), k))
            # Diagonale bas-gauche
            if i + 1 < N and j - 1 >= 0:
                springs.append(Spring(current_idx, index(i+1, j-1), sqrt(2), k))

    return masses, springs

def clone_system(masses, springs):
    # Clone des masses (copie profonde des objets)
    new_masses = [Mass(np.copy(m.position), m.mass) for m in masses]
    for i, m in enumerate(masses):
        new_masses[i].velocity = np.copy(m.velocity)
        new_masses[i].prev_position = np.copy(m.prev_position)
        new_masses[i].v_half = np.copy(m.v_half)
        new_masses[i].fixed = m.fixed
    # Clone des ressorts (identique)
    new_springs = copy.deepcopy(springs)
    return new_masses, new_springs

def serialize_system(masses, springs):
    return (
        [mass_to_dict(m) for m in masses],
        [spring_to_dict(s) for s in springs]
    )

def mass_to_dict(m):
    return {
        'position': m.position.tolist(),
        'velocity': m.velocity.tolist(),
        'prev_position': m.prev_position.tolist(),
        'v_half': m.v_half.tolist(),
        'mass': m.mass,
        'fixed': m.fixed
    }

def spring_to_dict(s):
    return {
        'mass1_idx': s.mass1_idx,
        'mass2_idx': s.mass2_idx,
        'rest_length': s.rest_length,
        'stiffness': s.stiffness
    }

def deserialize_system(mass_dicts, spring_dicts):
    masses = []
    for d in mass_dicts:
        m = Mass(d['position'], d['mass'])
        m.velocity = np.array(d['velocity'])
        m.prev_position = np.array(d['prev_position'])
        m.v_half = np.array(d['v_half'])
        m.fixed = d['fixed']
        masses.append(m)
    springs = [Spring(**s) for s in spring_dicts]
    return masses, springs



def compute_forces(masses, springs, damping, gravity):

    # Réinitialise les forces
    for mass in masses:
        mass.force[:] = 0.0

    # Applique les forces des ressorts
    for spring in springs:
        m1 = masses[spring.mass1_idx]
        m2 = masses[spring.mass2_idx]

        delta = m2.position - m1.position
        distance = np.linalg.norm(delta)
        if distance == 0:
            continue

        direction = delta / distance
        force_magnitude = spring.stiffness * (distance - spring.rest_length)
        force = force_magnitude * direction

        m1.force += force
        m2.force -= force

    # Gravité + amortissement visqueux
    # for mass in masses:
    #     if not mass.fixed:
    #         mass.force += np.array([0.0, gravity * mass.mass])     # gravité
    #         mass.force -= damping * mass.velocity                  # amortissement


def euler_explicit(masses, springs, dt, damping, gravity):
    compute_forces(masses, springs, damping, gravity)
    for mass in masses:
        if mass.fixed: continue
        acc = mass.force / mass.mass
        mass.velocity += acc * dt
        mass.position += mass.velocity * dt

def euler_implicit(masses, springs, dt, damping, gravity):
    compute_forces(masses, springs, damping, gravity)
    for mass in masses:
        if mass.fixed: continue
        acc = mass.force / mass.mass
        mass.velocity += acc * dt
        mass.position += mass.velocity * dt

def verlet_pos(masses, springs, dt, damping, gravity):
    compute_forces(masses, springs, damping, gravity)
    for mass in masses:
        if mass.fixed: continue
        acc = mass.force / mass.mass
        new_pos = 2 * mass.position - mass.prev_position + acc * dt**2
        mass.prev_position = np.copy(mass.position)
        mass.position = new_pos

def verlet_vel(masses, springs, dt, damping, gravity):
    # Première moitié
    compute_forces(masses, springs, damping, gravity)
    for mass in masses:
        if mass.fixed: continue
        acc = mass.force / mass.mass
        mass.position += mass.velocity * dt + 0.5 * acc * dt**2
        mass._acc_prev = acc  # temporairement stockée

    # Nouvelle accélération
    compute_forces(masses, springs, damping, gravity)
    for mass in masses:
        if mass.fixed: continue
        acc = mass.force / mass.mass
        mass.velocity += 0.5 * (mass._acc_prev + acc) * dt

def leapfrog(masses, springs, dt, damping, gravity):
    compute_forces(masses, springs, damping, gravity)
    for mass in masses:
        if mass.fixed: continue
        # Position
        mass.position += mass.v_half * dt

    compute_forces(masses, springs, damping, gravity)
    for mass in masses:
        if mass.fixed: continue
        acc = mass.force / mass.mass
        mass.v_half += acc * dt
        mass.velocity = mass.v_half - 0.5 * acc * dt  # pour lecture éventuelle

def rk4(masses, springs, dt, damping, gravity):
    positions0 = [np.copy(m.position) for m in masses]
    velocities0 = [np.copy(m.velocity) for m in masses]

    def compute_acc(pos_list):
        # remplace les positions temporaires
        for m, pos in zip(masses, pos_list):
            m.position = pos
        compute_forces(masses, springs, damping, gravity)
        return [m.force / m.mass for m in masses]

    a1 = compute_acc(positions0)
    v1 = velocities0

    positions2 = [p + 0.5*dt*v for p, v in zip(positions0, v1)]
    velocities2 = [v + 0.5*dt*a for v, a in zip(velocities0, a1)]
    a2 = compute_acc(positions2)

    positions3 = [p + 0.5*dt*v for p, v in zip(positions0, velocities2)]
    velocities3 = [v + 0.5*dt*a for v, a in zip(velocities0, a2)]
    a3 = compute_acc(positions3)

    positions4 = [p + dt*v for p, v in zip(positions0, velocities3)]
    velocities4 = [v + dt*a for v, a in zip(velocities0, a3)]
    a4 = compute_acc(positions4)

    for i, mass in enumerate(masses):
        if mass.fixed: continue
        mass.position = positions0[i] + (dt/6) * (v1[i] + 2*velocities2[i] + 2*velocities3[i] + velocities4[i])
        mass.velocity = velocities0[i] + (dt/6) * (a1[i] + 2*a2[i] + 2*a3[i] + a4[i])

def compute_energy(masses, springs, gravity=-0.1):
    E_kinetic = 0.0
    E_spring = 0.0
    E_gravity = 0.0

    for mass in masses:
        if not mass.fixed:
            v2 = np.dot(mass.velocity, mass.velocity)
            E_kinetic += 0.5 * mass.mass * v2
            E_gravity += mass.mass * gravity * mass.position[1]  # potentiel gravité

    for spring in springs:
        m1 = masses[spring.mass1_idx]
        m2 = masses[spring.mass2_idx]
        delta = m2.position - m1.position
        length = np.linalg.norm(delta)
        delta_length = length - spring.rest_length
        E_spring += 0.5 * spring.stiffness * delta_length**2

    return E_kinetic + E_spring + E_gravity

def process_system(data):
    index, mass_dicts, spring_dicts, integrator_name, dt, damping, gravity = data
    masses, springs = deserialize_system(mass_dicts, spring_dicts)

    # Récupération de la fonction à partir de son nom
    integrator_func = {
        "euler_explicit": euler_explicit,
        "euler_implicit": euler_implicit,
        "verlet_pos": verlet_pos,
        "verlet_vel": verlet_vel,
        "rk4": rk4,
        "leapfrog": leapfrog
    }[integrator_name]

    integrator_func(masses, springs, dt, damping, gravity)
    energy = compute_energy(masses, springs, gravity)
    return index, serialize_system(masses, springs), energy

def animate_mass_spring_systems(systems, titles, N, M, integrators, integrator_names, dt=0.01, damping=1.0, gravity=0.0, steps=200):
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    spring_lines = []
    mass_dots = []

    for ax, (masses, springs), title in zip(axes, systems, titles):
        ax.set_title(title)
        ax.set_aspect('equal')
        lines = []
        for spring in springs:
            p1 = masses[spring.mass1_idx].position
            p2 = masses[spring.mass2_idx].position
            line, = ax.plot([p1[0], p2[0]], [p1[1], p2[1]], 'k-', lw=1)
            lines.append(line)
        spring_lines.append(lines)

        dots, = ax.plot(
            [m.position[0] for m in masses],
            [m.position[1] for m in masses],
            'ro', markersize=3
        )
        ax.set_xlim([-1, N + 1])
        ax.set_ylim([-1, M + 1])
        mass_dots.append(dots)

    # Fenêtre pour la courbe d’énergie
    fig_energy, ax_energy = plt.subplots(figsize=(10, 4))
    ax_energy.set_title("Évolution de l'énergie totale")
    ax_energy.set_xlabel("Temps")
    ax_energy.set_ylabel("Énergie")
    ax_energy.set_xlim(0, steps)

    energy_lines = []
    energy_histories = [[] for _ in systems]
    energy_refs = []

    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']

    for i, (masses, springs) in enumerate(systems):
        energy = compute_energy(masses, springs, gravity)
        energy_refs.append(energy)
        line, = ax_energy.plot([], [], label=titles[i], color=colors[i % len(colors)])
        ax_energy.axhline(y=energy, color=colors[i % len(colors)], linestyle="--", alpha=0.3)
        energy_lines.append(line)

    ax_energy.legend()

    def update(frame):
        task_data = []
        for i, (masses, springs) in enumerate(systems):
            mass_dicts, spring_dicts = serialize_system(masses, springs)
            task_data.append((i, mass_dicts, spring_dicts, integrator_names[i], dt, damping, gravity))

        with mp.Pool() as pool:
            results = pool.map(process_system, task_data)

        for i, (mass_dicts, spring_dicts), energy in results:
            masses, springs = deserialize_system(mass_dicts, spring_dicts)
            systems[i] = (masses, springs)

            # Mise à jour visualisation des ressorts
            for j, spring in enumerate(springs):
                p1 = masses[spring.mass1_idx].position
                p2 = masses[spring.mass2_idx].position
                spring_lines[i][j].set_data([p1[0], p2[0]], [p1[1], p2[1]])

            # Mise à jour des masses
            mass_dots[i].set_data(
                [m.position[0] for m in masses],
                [m.position[1] for m in masses]
            )
            axes[i].set_title(f"{titles[i]}\nÉnergie: {energy:.3f}")

            # Mise à jour de l'historique énergétique
            energy_histories[i].append(energy)
            energy_lines[i].set_data(range(len(energy_histories[i])), energy_histories[i])

        # Mise à jour des limites dynamiques du graphe d'énergie
        ax_energy.set_xlim(0, max(steps, max(len(hist) for hist in energy_histories)))
        max_energy = max((max(hist) for hist in energy_histories if hist), default=1)
        ax_energy.set_ylim(0, max_energy * 1.1)

        return sum(spring_lines, []) + mass_dots + energy_lines
    # Animation principale
    ani1 = animation.FuncAnimation(
        fig, update, frames=steps, interval=30, blit=False, repeat=False
    )

    # Fonction pour mise à jour de la figure énergétique
    def update_energy(frame):
        for i, line in enumerate(energy_lines):
            line.set_data(range(len(energy_histories[i])), energy_histories[i])
        max_energy = max((max(hist) for hist in energy_histories if hist), default=1)
        min_energy = min((min(hist) for hist in energy_histories if hist), default=1)
        ax_energy.set_xlim(0, max(100, max(len(hist) for hist in energy_histories)))
        ax_energy.set_ylim(min_energy, max_energy)
        print(min_energy, max_energy)
        return energy_lines

    ani2 = animation.FuncAnimation(
        fig_energy, update_energy, frames=steps, interval=30, blit=False, repeat=False
    )

    plt.suptitle("Évolution dynamique des méthodes d'intégration", fontsize=16)
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.show()



if __name__ == "__main__":
    # Création initiale
    N, M = 5,5
    k = 5.0
    dt = 0.5
    damping = 1.0
    gravity = -0.1 #-0.05
    masses_base, springs_base = create_mass_spring_grid(N, M, k)

    # Clonage
    displayed = [True, True, True, True, True, True]
    methods = ["Euler explicite", "Euler implicite", "Verlet position", "Verlet vitesse", "RK4", "Leapfrog"]
    integrator_names = ["euler_explicit", "euler_implicit", "verlet_pos", "verlet_vel", "rk4", "leapfrog"]
    integrators = [euler_explicit, euler_implicit, verlet_pos, verlet_vel, rk4, leapfrog]
    systems = []

    for integrator in integrators:
        m_clone, s_clone = clone_system(masses_base, springs_base)
        integrator(m_clone, s_clone, dt, damping, gravity)
        systems.append((m_clone, s_clone))

    # Affichage des 6 versions
    animate_mass_spring_systems(
        [system for i, system in enumerate(systems) if displayed[i]],
        [method for i, method in enumerate(methods) if displayed[i]],
        N, M,
        integrators,
        integrator_names,
        dt, damping, gravity
    )

