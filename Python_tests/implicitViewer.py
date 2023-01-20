import math

import matplotlib.pyplot as plt
import numpy as np

def hessian(x):
    """
    Calculate the hessian matrix with finite differences
    Parameters:
       - x : ndarray
    Returns:
       an array of shape (x.dim, x.ndim) + x.shape
       where the array[i, j, ...] corresponds to the second derivative x_ij
    """
    x_grad = np.gradient(x)
    hessian = np.empty((x.ndim, x.ndim) + x.shape, dtype=x.dtype)
    for k, grad_k in enumerate(x_grad):
        # iterate over dimensions
        # apply gradient again to every component of the first derivative.
        tmp_grad = np.gradient(grad_k)
        for l, grad_kl in enumerate(tmp_grad):
            hessian[k, l, :, :] = grad_kl
    return hessian

def wyvill(x):
    return 1 - ((1 - x)**2)**3
delta = 0.0025
xrange = np.arange(0, 1, delta)
yrange = np.arange(0, 1, delta)
X, Y = np.meshgrid(xrange, yrange)

Y = 1 - Y

sqrt2 = math.sqrt(2)
n = 2
# F = (X**n + Y**n)**(1/n) / sqrt2
def h(x, y):
    w = np.clip(1 - ((x + y - 1.0) / (2.0 * y-1.0))**(1.0/(1.0-y)), 0, 1)
    return w

def F(x, y):
    # return np.maximum(x, y)
    # return (x**n + y**n)**(1/n)
    # return (y <= 0.5) * ((x <= 0.5) * np.maximum(x, y) + (x > 0.5) * (1 - np.maximum(x / 2, y / 2))) + (y > 0.5) * y #  Find the good replacement function
    # return ((np.cos(x * 10)*0.5 + 0.8)**n + (np.sin(y*10)*0.5 + 0.5)**n)**(1/n)
    return ((x + y) > 1.0) * ((0.5 + (np.maximum(x, y) - 0.5)) * h(np.minimum(x, y), np.maximum(x, y))) + ((x + y) <= 1.0) * np.maximum(x, y)
    # return (np.maximum(x, y) <= 0.5) * np.maximum(x, y) + (np.maximum(x, y) > 0.5) * (np.abs(x - y) + 0.5)

def grad_(x, y):
    return np.gradient(F(x, y))


plot = F(X, Y)
grad = grad_(X, Y)  # np.gradient(plot)

_dxF = abs(grad[0])
_dyF = abs(grad[1])

dxF = _dxF / np.sqrt(_dxF**2 + _dyF**2)
dyF = _dyF / np.sqrt(_dxF**2 + _dyF**2)

B = (X*X)/(X*X + Y*Y)

# A = dxF[200, 200]

# plt.imshow(dxF, cmap="seismic", vmin=0, vmax=1)
# plt.imshow(dxF / (dxF + dyF), cmap="seismic", vmin=0, vmax=1)
# plt.imshow(B - (dyF / (dxF + dyF)), cmap="seismic")
# plt.imshow(B, cmap="seismic", vmin=0, vmax=1)
# plt.colorbar()

plt.imshow(plot, cmap="seismic", vmin=0, vmax=1)
plt.contour(plot, [0.5])
# plt.clabel(plt.contour(plot, [0.5, 0.6, 0.7, 0.8, 0.9]))
# plt.title("Material for $f=(x^n + y^n)^{1/n}$ with $n=" + str(n) + "$")
plt.title("Blend: $(x^n + y^n)^{1/n}$ with $n=" + str(n) + "$")
# plt.title("Blend = $max(x, y)$")
# plt.title("Materials")
# plt.title("Blend: Simple contact")
#
# def f(X, Y):
#     Ma = np.clip(1 - (np.sqrt((0.25 - X) ** 2 + (0.5 - Y) ** 2) - 0.2) - 0.5, 0, 1)
#     Mb = np.clip(1 - (np.sqrt((0.75 - X) ** 2 + (0.5 - Y) ** 2) - 0.2) - 0.5, 0, 1)
#     # return (Ma + Mb)/2
#     g = grad_(X, Y)
#     _dxF = abs(g[0])
#     _dyF = abs(g[1])
#
#     dxF = _dxF / np.sqrt(_dxF**2 + _dyF**2)
#     dyF = _dyF / np.sqrt(_dxF**2 + _dyF**2)
#
#     img = np.clip(F(Ma, Mb) - 0.5, 0, 1)
#
#     return (dxF) * (1 - img) * (F(Ma, Mb) > 0.5)
#     # return ((dyF) * (img) + (dxF) * (1 - img) ) * (F(Ma, Mb) > 0.5)
#
# M = f(X, Y)
# plt.imshow(M, cmap="seismic", vmin=0, vmax=1)
# # plt.contour(M, [0.5]) #np.arange(0, 1, 0.01))
plt.figaspect(1)
plt.xticks([])
plt.yticks([])
plt.show()


