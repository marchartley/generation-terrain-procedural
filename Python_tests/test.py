import math
import numpy as np
import matplotlib.pyplot as plt


T_implicit = np.zeros((100, 100, 100))
T_voxels = np.zeros((100, 100, 100))
T_heightmap = np.zeros((100, 100))

Q = 1.0
C = np.array([50, 50, 50])
R = 10

for x in range(T_implicit.shape[0]):
    for y in range(T_implicit.shape[1]):
        p = np.array([x, y, 50])
        d = np.linalg.norm(p - C)
        if d < R:
            T_heightmap[x, y] = (Q / (R ** 3 * math.pi * 2 / 3)) * math.sqrt(R ** 2 - d ** 2)
        for z in range(T_implicit.shape[2]):
            p = np.array([x, y, z])
            d = np.linalg.norm(p - C)
            if d < R:
                T_implicit[x, y, z] = (3.0 / (math.pi * R**3)) * Q * (1 - (d/R))
                T_voxels[x, y, z] = (1 - (d/R))

T_voxels = T_voxels * Q / np.sum(T_voxels)
print(np.sum(T_implicit))
print(np.sum(T_voxels))
print(np.sum(T_heightmap))


plt.imshow(T_implicit[:, :, 50])
plt.show()
plt.imshow(T_voxels[:, :, 50])
plt.show()
plt.imshow(T_heightmap)
plt.show()


for x in range(T_implicit.shape[0]):
    for y in range(T_implicit.shape[1]):
        for z in range(T_implicit.shape[2]):
            T_implicit[x, y, 0] += T_implicit[x, y, z]
print(np.sum(T_implicit[:, :, 0]))
plt.imshow(T_implicit[:, :, 0])
plt.show()