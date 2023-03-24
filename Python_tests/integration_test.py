import numpy as np
import matplotlib.pyplot as plt

m = 1.0  # mass of particle
k = 3.0  # spring constant

# Define acceleration function
def acceleration(x):
    return -k * x / m

# Example code for Euler integration

def euler_integration(x, v, a, dt):
    # Update velocity
    v += a * dt

    # Update position
    x += v * dt

    return x, v


# Example code for Verlet integration

def verlet_integration(x, x_prev, a, dt):
    # Update position
    x_next = 2 * x - x_prev + a * dt ** 2

    return x_next, x


# Example code for Leapfrog Verlet integration

def leapfrog_verlet_integration(x, v, a, dt):
    # Update velocity halfway
    v_half = v + a * dt / 2

    # Update position
    x += v_half * dt

    # Update acceleration using new position
    a_new = acceleration(x)

    # Update velocity using new acceleration
    v_new = v_half + a_new * dt / 2

    return x, v_new


# Define function for exact position
def exact_position(x0, v0, t):
    omega = np.sqrt(k / m) # angular frequency
    return x0 * np.cos(omega * t) + v0 / omega * np.sin(omega * t)



def main():
    # Define parameters
    x0 = 1.0  # initial position
    v0 = 0.0  # initial velocity
    dt = 0.10  # time step
    t_max = 100.0  # maximum simulation time

    # Define initial conditions
    x_euler = x0
    v_euler = v0
    x_verlet = x0
    x_prev_verlet = x0 - v0 * dt
    x_leapfrog = x0
    v_leapfrog = v0

    # Define arrays to store positions and times
    t_array = np.arange(0, t_max, dt)
    x_array_euler = np.zeros_like(t_array)
    x_array_verlet = np.zeros_like(t_array)
    x_array_leapfrog = np.zeros_like(t_array)

    # Loop over time steps and update positions using each integration method
    for i in range(len(t_array)):
        # Euler integration
        x_euler, v_euler = euler_integration(x_euler, v_euler, acceleration(x_euler), dt)
        x_array_euler[i] = x_euler

        # Verlet integration
        x_verlet, x_prev_verlet = verlet_integration(x_verlet, x_prev_verlet, acceleration(x_verlet), dt)
        x_array_verlet[i] = x_verlet

        # Leapfrog Verlet integration
        x_leapfrog, v_leapfrog = leapfrog_verlet_integration(x_leapfrog, v_leapfrog, acceleration(x_leapfrog), dt)
        x_array_leapfrog[i] = x_leapfrog

    # Calculate exact positions at each time step
    x_exact = exact_position(x0, v0, t_array)

    # Plot exact position and compare to numeric integrations
    plt.plot(t_array, x_exact, label='Exact')
    # plt.plot(t_array, x_array_euler, label='Euler')
    # plt.plot(t_array, x_array_verlet, label='Verlet')
    plt.plot(t_array, x_array_leapfrog, label='Leapfrog Verlet')
    plt.legend()
    plt.xlabel('Time')
    plt.ylabel('Position')
    plt.show()


if __name__ == "__main__":
    main()