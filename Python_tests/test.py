import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from PIL import Image

# Load the image and convert it to grayscale
image_path = "lena.jpg"
image = Image.open(image_path).convert("L").resize((1000, 1000))  # "L" mode converts to grayscale
image = image.crop((0, 500, 1000, 800))
image_array = np.array(image)

# Get dimensions of the image
height, width = image_array.shape

# Create X, Y coordinate arrays
x = np.arange(0, width, 1)
y = np.arange(0, height, 1)
X, Y = np.meshgrid(x, y)

# Create Z array with intensity values from the image
Z = image_array

# Plot the 3D surface
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot surface
surf = ax.plot_surface(X, Y, Z, cmap=cm.gray, linewidth=1, antialiased=True)

# Customize the z axis.
ax.set_zlim(0, 516)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
