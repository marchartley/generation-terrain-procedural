import numpy as np
import matplotlib.pyplot as plt

class point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

img_file = "C:/Users/march/Desktop/MNT_MAY100m_HOMONIM_WGS84_NM_ZNEG.glz"

x = []
y = []
z = []
zs = []
with open(img_file) as file:
    for line in file:
        s_line = line.split(" ")
        orig_x, orig_y, orig_z = s_line
        _x, _y, _z = round(float(orig_x) * 1000), round(float(orig_y) * 1000), float(orig_z)
        x.append(_x)
        y.append(_y)
        z.append(_z)

minX, maxX = min(x), max(x)
minY, maxY = min(y), max(y)
minZ, maxZ = min(z), max(z)

array = np.zeros((maxX - minX +1, maxY - minY +1))

with open(img_file) as file:
    for line in file:
        s_line = line.split(" ")
        orig_x, orig_y, orig_z = s_line
        _x, _y, _z = round(float(orig_x) * 1000), round(float(orig_y) * 1000), float(orig_z)
        array[_x - minX, (maxY - minY) - (_y - minY)] = _z

plt.imshow(array.swapaxes(0, 1))
plt.show()
