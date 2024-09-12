import numpy as np
import matplotlib.pyplot as plt
from scipy.special import k0

# Define the grid size and parameters
grid_size = 50
x0, y0 = 25, 25  # Source location
D = 1.0  # Diffusion coefficient
lambda_decay = 0.01  # Decay rate
S0 = 100  # Emission strength
u, v = 0, 0  # Advection velocities in x and y directions

# Create a meshgrid for the coordinates
x = np.linspace(0, grid_size-1, grid_size)
y = np.linspace(0, grid_size-1, grid_size)
X, Y = np.meshgrid(x, y)

# Calculate distances from the source point
distances = np.sqrt((X - x0)**2 + (Y - y0)**2)
distances[distances == 0] = 1e-10  # Avoid division by zero at the source

# Calculate concentration using the modified Bessel function
concentration = S0 * k0(np.sqrt(lambda_decay / D) * distances) / (2 * np.pi * D)

# Adjust for advection using central difference for the gradient
dx = dy = 1  # Grid spacing
grad_c_x = np.gradient(concentration, dx, axis=1)
grad_c_y = np.gradient(concentration, dy, axis=0)

# Compute advection term
advection_effect = -u * grad_c_x - v * grad_c_y
concentration += advection_effect  # Add the advection effect to the concentration

# Plot the concentration field
plt.figure(figsize=(10, 8))
plt.imshow(concentration, cmap='hot', origin='lower')
plt.colorbar(label='Concentration')
plt.title('Concentration at Equilibrium with Advection')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.show()




# import numpy as np
# from scipy.optimize import minimize
# from scipy.special import k0
# import matplotlib.pyplot as plt
#
# # Target values
# Q = 99000  # Desired integral
# target_c0 = 10.0  # Desired concentration at the source
#
# # Define the grid size and parameters
# grid_size = 500
# x0, y0 = 250, 250  # Source location
# S0 = 100  # Emission strength
# epsilon = 5  # Small radius around the source for regularization
#
# # Create a meshgrid for the coordinates
# x = np.linspace(0, grid_size - 1, grid_size)
# y = np.linspace(0, grid_size - 1, grid_size)
# X, Y = np.meshgrid(x, y)
# # Calculate distances from the source point
# distances = np.sqrt((X - x0)**2 + (Y - y0)**2)
# # Regularization: Adjust values close to the source
# within_epsilon = distances < epsilon
#
#
#
# def objective_function(params):
#     D, lambda_decay = params
#
#     average_value = k0(np.sqrt(lambda_decay / D) * epsilon) / (2 * np.pi * D) * S0
#     concentration = np.where(within_epsilon, average_value, S0 * k0(np.sqrt(lambda_decay / D) * distances) / (2 * np.pi * D))
#
#     # Compute the integral of the concentration
#     total_concentration = np.sum(concentration)
#     return (total_concentration - Q)**2 + (concentration[x0, y0] - target_c0)**2
#
#
# # Initial guesses for D and lambda_decay
# initial_guess = [1.0, 0.01]
#
# result = minimize(objective_function, initial_guess, method='Nelder-Mead')
#
# optimized_D, optimized_lambda_decay = -1, -1
#
# # if result.success:
# optimized_D, optimized_lambda_decay = result.x
# print(f"Optimized D: {optimized_D}")
# print(f"Optimized Lambda Decay: {optimized_lambda_decay}")
# D = optimized_D
# lambda_decay = optimized_lambda_decay
# average_value = k0(np.sqrt(lambda_decay / D) * epsilon) / (2 * np.pi * D) * S0
# concentration = np.where(within_epsilon, average_value, S0 * k0(np.sqrt(lambda_decay / D) * distances) / (2 * np.pi * D))
#
# # Compute the integral of the concentration
# total_concentration = np.sum(concentration)
# print(f"Total concentration in the field: {total_concentration}")
# # Plot the concentration field
# plt.figure(figsize=(10, 8))
# plt.imshow(concentration, cmap='hot', origin='lower')
# plt.colorbar(label='Concentration')
# plt.title('Concentration at Equilibrium')
# plt.xlabel('X Coordinate')
# plt.ylabel('Y Coordinate')
# plt.show()
# # else:
# #     print("Optimization failed.")
# #     print(result.x)
#
#
#
#
# # import numpy as np
# # import matplotlib.pyplot as plt
# # from scipy.special import k0
# #
# # # Define the grid size and parameters
# # grid_size = 500
# # x0, y0 = 250, 250  # Source location
# # D = 1.0  # Diffusion coefficient
# # lambda_decay = 0.01  # Decay rate
# # S0 = 100  # Emission strength
# # epsilon = 5  # Small radius around the source for regularization
# #
# # # Create a meshgrid for the coordinates
# # x = np.linspace(0, grid_size - 1, grid_size)
# # y = np.linspace(0, grid_size - 1, grid_size)
# # X, Y = np.meshgrid(x, y)
# #
# # # Calculate distances from the source point
# # distances = np.sqrt((X - x0)**2 + (Y - y0)**2)
# #
# # # Regularization: Adjust values close to the source
# # within_epsilon = distances < epsilon
# # average_value = k0(np.sqrt(lambda_decay / D) * epsilon) / (2 * np.pi * D) * S0
# # concentration = np.where(within_epsilon, average_value, S0 * k0(np.sqrt(lambda_decay / D) * distances) / (2 * np.pi * D))
# #
# # # Plot the concentration field
# # plt.figure(figsize=(10, 8))
# # plt.imshow(concentration, cmap='hot', origin='lower')
# # plt.colorbar(label='Concentration')
# # plt.title('Concentration at Equilibrium')
# # plt.xlabel('X Coordinate')
# # plt.ylabel('Y Coordinate')
# # plt.show()
# #
# # # Compute the integral of the concentration
# # total_concentration = np.sum(concentration)
# # print(f"Total concentration in the field: {total_concentration}")
