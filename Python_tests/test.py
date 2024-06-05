import numpy as np
import matplotlib.pyplot as plt

def sediment_concentration(x, y, x0, y0, ux, uy, D, t):
    """
    Calculate the concentration of sediments at point (x, y) after time t.
    """
    return (1 / (4 * np.pi * D * t)) * np.exp(-((x - x0 - ux * t)**2 + (y - y0 - uy * t)**2) / (4 * D * t))

# Parameters
x0, y0 = 0.0, 0.0   # Initial point
ux, uy = 1.0, 0.5   # Velocity field
D = 1.0             # Diffusion coefficient
t = 1.0             # Time

# Grid definition
x = np.linspace(-10, 10, 400)
y = np.linspace(-10, 10, 400)
X, Y = np.meshgrid(x, y)

# Compute concentration
Z = sediment_concentration(X, Y, x0, y0, ux, uy, D, t)

# Plotting
plt.figure(figsize=(8, 6))
plt.contourf(X, Y, Z, levels=50, cmap='viridis')
plt.colorbar(label='Concentration')
plt.title(f'Sediment Concentration at t={t}')
plt.xlabel('x')
plt.ylabel('y')
plt.scatter([x0 + ux * t], [y0 + uy * t], color='red')  # Indicate the advected source position
plt.annotate('Source Position', (x0 + ux * t, y0 + uy * t), textcoords="offset points", xytext=(5,-10), ha='center', color='red')
plt.grid(True)
plt.show()
