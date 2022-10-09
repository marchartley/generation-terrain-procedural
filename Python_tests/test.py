import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = plt.plot([], [], 'r-')

def init():
    # ax.set_xlim(0, 2*np.pi)
    # ax.set_ylim(-1, 1)
    ax.set_aspect('equal')
    return ln,


def gaussian(x, mu, sig, wrap = False):
    if wrap:
        return gaussian(x, mu -1, sig) + gaussian(x, mu, sig) + gaussian(x, mu + 1, sig)
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

mu = 0.0
sigma = 0.1
def update(_mu):
    if _mu > .5:
        ani.event_source.stop()
        plt.close(fig)
        return
    global mu, sigma
    mu = _mu
    vals = np.linspace(0, 1, 128)
    x = np.cos(np.linspace(0, 2*np.pi, 128))
    y = np.sin(np.linspace(0, 2*np.pi, 128))
    r = gaussian(vals, mu - int(mu), sigma, True) + gaussian(vals, mu - int(mu) - .5, sigma, True)
    # plt.plot(vals, r); return
    # y = r
    # x = vals

    x *= r
    y *= r

    ln.set_data(x, y)
    ax.set_xlim(min(x), max(x))
    ax.set_ylim(min(y), max(y))
    # mu += .1
    # sigma += .01

    return ln,

# update(None)
ani = FuncAnimation(fig, update, frames = np.linspace(0, 1, 50),
                    init_func=init, blit=False, interval=50)
plt.show()