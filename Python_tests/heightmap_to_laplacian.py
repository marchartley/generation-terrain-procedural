import PIL
from PIL import Image
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def gradient(img: np.ndarray) -> np.ndarray:
    # res = np.zeros((img.shape[0], img.shape[1], 2))
    # for y in range(1, img.shape[0] - 1):
    #     for x in range(1, img.shape[1] - 1):
    #         grad = np.array([
    #             (img[y + 1, x] - img[y - 1, x]) * .5,
    #             (img[y, x + 1] - img[y, x - 1]) * .5
    #         ])
    #         res[y, x] = grad
    # return res
    [yGrad, xGrad] = np.gradient(img)
    grad = np.zeros((xGrad.shape[0], xGrad.shape[1], 2))
    grad[:, :, 0] = yGrad
    grad[:, :, 1] = xGrad
    return grad

def diverence(img: np.ndarray) -> np.ndarray:
    res = np.zeros((img.shape[0], img.shape[1]))
    for y in range(1, img.shape[0] - 1):
        for x in range(1, img.shape[1] - 1):
            div = ((img[y + 1, x, 1] - img[y - 1, x, 1]) + (img[y, x + 1, 0] - img[y, x - 1, 0])) * .5
            res[y, x] = div
    return res

def main():
    # heightmap_filename = r"diffusion_curve_result.png"
    heightmap_filename = r"C:\codes\Qt\Shared_folder_IRIT_Seafile\generation-terrain-procedural\saved_maps\heightmaps\___new_one_slope.png"
    heightmap = np.array(Image.open(heightmap_filename).resize((100, 100), resample=Image.LINEAR))[:,:,0].astype(float)
    # plt.show()
    grads = gradient(heightmap)
    grads /= np.max(grads)
    lap = diverence(grads)
    i_grads = np.zeros((grads.shape[0], grads.shape[1], 3))
    i_grads[:, :, :2] = abs(grads)
    i_grads /= np.max(i_grads)
    absMax = max(abs(np.min(lap)), abs(np.max(lap)))
    lap = 1 + 0.5 * ((lap - absMax)/ absMax) # / (np.max(lap) - np.min(lap)) #255.0 - ((lap * 255.0/2.0) + 128.0)
    fig, ax = plt.subplots(1, 3)
    ax[0].imshow(heightmap, cmap="gray")
    ax[1].imshow(i_grads)
    asRGB = np.zeros((lap.shape[0], lap.shape[1], 3))
    asRGB[:, :, 0] = asRGB[:, :, 1] = asRGB[:, :, 2] = lap
    ax[2].imshow(asRGB)
    plt.imsave("lap.png", asRGB)
    plt.show()


if __name__ == "__main__":
    main()
