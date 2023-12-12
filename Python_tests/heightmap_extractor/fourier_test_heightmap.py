import math
import random
import time

import PIL.Image
from scipy.fft import fft, ifft
import numpy as np
from PIL import Image
from math import floor
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import noise.perlin as perlin
from scipy import signal

def clamp(x, _min, _max):
    return _min if x < _min else _max if x > _max else x

def mutate_heights(y : np.ndarray, returned_array: np.ndarray, axis: int):

    mu, sigma = 1., 0.1
    phase_mu, phase_sigma = 0.0, 0.1
    noise = perlin.SimplexNoise()
    offset1 = 10
    offset2 = 100

    gradients = [0] + [y[i] - y[i - 1] for i in range(1 , len(y))]
    grad2D = np.array([gradients for _ in range(len(gradients))]) * 0
    # grad2D = (np.rot90(grad2D, 0) + np.rot90(grad2D, 1)+ np.rot90(grad2D, 2)+ np.rot90(grad2D, 3)) * 10

    widths = [i for i in range(1, 1 + len(y))]
    cwtmatr = signal.cwt(y, signal.wavelets.ricker, widths)
    # cwtmatr = signal.cwt(cwtmatr[:, 0], signal.wavelets.ricker, widths)
    final = (np.rot90(cwtmatr, 0)  + np.rot90(cwtmatr, 1)  + np.rot90(cwtmatr, 2)  + np.rot90(cwtmatr, 3)) + grad2D

    scale = .01
    nosy = np.array([[noise.noise2(x * scale, y * scale) for x in range(len(gradients))] for y in range(len(gradients))])

    return final * nosy

    for i in range(len(y)):
        i_x = (i/len(y)) * 2.0 + offset1
        i_y = (i/len(y)) * 2.0 + offset2
        fourier = fft(np.roll(y, 0))
        phase_shift = 0.0 + (noise.noise2(i_x, i_y) * 0.1) #random.gauss(phase_mu, phase_sigma)
        ampli_shift = 1.0 + (noise.noise2(i_x, i_y) * 1.1) #random.gauss(mu, sigma)
        # t_fourier = [f * complex(math.cos(phase_shift), math.sin(phase_shift)) * ampli_shift * complex(noise.noise2(i_f, i_f), noise.noise2(i_f, i_f)) for f in fourier]
        t_fourier = [f * complex(math.cos(phase_shift), math.sin(phase_shift)) * ampli_shift for n, f in enumerate(fourier)]
        # t_fourier = [f for f in fourier]
        # t_fourier = [f * complex(abs(noise.noise2((i + 10000*axis)/100, (100 + i + 100000*axis))/100), 1.) for f in fourier]
        # t_fourier = [f * complex(1., random.gauss(mu, sigma)) for f in fourier]
        # t_fourier = [f * complex(random.gauss(mu, sigma), random.gauss(mu, sigma)) for f in fourier]

        y_modified = (ifft(t_fourier).real + 1.0) / 2
        if axis == 0:
            returned_array[i, :] = np.add(returned_array[i, :], y_modified)
        elif axis == 1:
            returned_array[:, i] = np.add(returned_array[:, i], y_modified)
        else:
            raise ValueError("Mauvais axe")
    return returned_array


def main():
    img: np.ndarray = np.asarray(Image.open("test_fourier_image.png").convert('L'))

    y = np.flip(img, 0).argmin(0)

    y = -((y / (.5 * img.shape[0])) - 1)

    img_size = img.shape[1]
    return_img = np.ones((img_size, img_size))

    return_img = mutate_heights(y, return_img, 0)
    # return_img = mutate_heights(y, return_img, 1)
    modified_img = PIL.Image.fromarray(return_img)
    # modified_img = modified_img.resize((50, 50))# .resize((100, 100))
    return_img = np.array(modified_img)
    normalized = (modified_img - return_img.min()) / (return_img.max() - return_img.min())
    modified_img = PIL.Image.fromarray(normalized * 255)
    modified_img.convert('RGB').save("result.jpg")

    # plt.imshow(return_img)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    X, Y = np.meshgrid(range(return_img.shape[0]), range(return_img.shape[1]))
    # ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1, 1, 0.5, 1]))
    ax.set_box_aspect(aspect=(1, 1, 0.2))
    ax.plot_surface(X, Y, return_img, cmap='viridis', edgecolor='none')
    plt.show()


if __name__ == "__main__":
    main()