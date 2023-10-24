import time
import numpy as np
import PIL.Image
import matplotlib.pyplot as plt

from numpy.fft  import fft2, ifft2

from coralize_my_island import *

def np_fftconvolve(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    return np.real(ifft2(fft2(A)*fft2(B, s=A.shape)))

def fft_convolve2d(x,y):
    """ 2D convolution, using FFT"""
    fr = fft2(x)
    fr2 = fft2(np.flipud(np.fliplr(y)), s=x.shape)
    m,n = fr.shape
    cc = np.real(ifft2(fr*fr2))
    # cc = np.roll(cc, -m//2+1,axis=0)
    # cc = np.roll(cc, -n//2+1,axis=1)
    cc = np.roll(cc, int(-y.shape[0]//2 + 1), axis=0)
    cc = np.roll(cc, int(-y.shape[1]//2 + 1), axis=1)
    return cc

def blur(img: np.ndarray, mask: np.ndarray, kernelSize: Union[Tuple[int, int], int] = (5, 5)) -> np.ndarray:
    ksize = kernelSize
    if isinstance(ksize, int):
        ksize = (ksize, ksize)
    kernel = np.ones(ksize) / (ksize[0] * ksize[1])
    blurredVersion = fft_convolve2d(img, kernel)
    blurredVersion[mask == 1] = img[mask == 1]
    return blurredVersion

def diffuseValues(initialMap: np.ndarray, mask: np.ndarray, kernelSize: Union[Tuple[int, int], int] = (5, 5),
                  epsilon: float = 1e-5, maxIterations: int = 10000,
                  verbose: Union[bool, int] = False) -> np.ndarray:
    if isinstance(verbose, bool):
        verbose = 1 if verbose else 0

    blurred = initialMap.copy()
    previousValues = initialMap.copy()
    nbIterations = 0
    difference = 10000  # Just a random high value
    for i in range(maxIterations):
        nbIterations = i + 1
        blurred = blur(blurred, mask, kernelSize=kernelSize)
        difference = np.mean(np.abs(blurred - previousValues))
        if verbose > 0 and i % verbose == 0:
            print(f"Iteration {nbIterations}: difference ~= 10^{int(math.log10(difference))}")
        if difference < epsilon:
            break
        previousValues = blurred.copy()

    if verbose > 0:
        print(f"Stopped after {nbIterations} iterations ", end = "")
        if difference < epsilon:
            print("(difference got very low)")
        else:
            print(f"(max iterations [{maxIterations}] reached)")
    return blurred

def main():
    initialImg = PIL.Image.open("skeleton_3.png")  # .resize((200, 200), resample=False)
    img = 1 - (np.array(initialImg)[:, :, 0] / 255)
    # Remove the color palette used in the drawing
    nbColorsUsed = len(np.unique(img)) - 1
    img[1:3, 1:(nbColorsUsed*2 + 1)] = 0

    # Raising to a power < 1 to keep values pretty close
    img = np.power(img, .1)
    mask = img.copy()
    mask[mask > 0] = 1

    startingTime = time.time()
    blurred = diffuseValues(img, mask, kernelSize = (3, 3), epsilon = 1e-8, verbose = False)
    print(f"Diffusion applied on {img.shape[1]}x{img.shape[0]} image in {round((time.time() - startingTime) * 1000)}ms")

    # One last blur to remove traces of the initial curves
    blurred = blur(blurred, np.zeros_like(mask), kernelSize = 3)

    plt.imsave("diffusion_curve_result.png", blurred, cmap="gray")

    plt.imshow(blurred)
    plt.show()

    plot3D(blurred)
    plt.show()


if __name__ == "__main__":
    main()
