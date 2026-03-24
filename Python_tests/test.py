import matplotlib.pyplot as plt
import numpy as np

m = 1.0
dt = 1.0
g = 1.0
T = 100.0
bounces = True
coef_rest = 0.5
point_gravity = False

P0 = np.array([20.0, 10.0])
V0 = np.array([10.0, 1.0])

center = np.array([20.0, 0.0])
def F(pos: np.ndarray) -> np.ndarray:
    if point_gravity:
        pc = (pos - center)
        r = np.linalg.norm(pc)
        G = (pc / r) * g
        return -G * m
    else:
        return np.array([0.0, -g])

eulerP = np.array([np.array(P0)])
eulerV = np.array(V0)
for i in range(int(T / dt)):
    p = np.array(eulerP[-1])
    v = np.array(eulerV)
    a = F(p) / m
    v += (a * dt)
    p += v * dt


    if bounces and p[1] < 0:
        v *= coef_rest
        v[1] *= -1.0
        p[1] = max(0, p[1])
    eulerP = np.concatenate([eulerP, np.array([p])], axis=0)
    eulerV = v


semiP = np.array([np.array(P0)])
semiV = np.array(V0)
for i in range(int(T / dt)):
    p = np.array(semiP[-1])
    v = np.array(semiV)
    a = F(p) / m
    p += v * dt
    v += (a * dt)

    if bounces and p[1] < 0:
        v *= coef_rest
        v[1] *= -1.0
        p[1] = max(0, p[1])
    semiP = np.concatenate([semiP, np.array([p])], axis=0)
    semiV = v


leapP = np.array([np.array(P0)])
leapV = np.array(V0)
for i in range(int(T / dt)):
    p = np.array(leapP[-1])
    v = np.array(leapV)
    a = F(p) / m
    p += v * dt * 0.5
    v += (a * dt)
    p += v * dt * 0.5

    if bounces and p[1] < 0:
        v *= coef_rest
        v[1] *= -1.0
        p[1] = max(0, p[1])
    leapP = np.concatenate([leapP, np.array([p])], axis=0)
    leapV = v


halfP = np.array([np.array(P0)])
halfV = np.array(V0)
for i in range(int(T / dt)):
    p = np.array(halfP[-1])
    v = np.array(halfV)
    a = F(p) / m
    v += (a * dt) * 0.5

    p += v * dt

    a = F(p) / m
    v += (a * dt) * 0.5

    if bounces and p[1] < 0:
        v *= coef_rest
        v[1] *= -1.0
        p[1] = max(0, p[1])
    halfP = np.concatenate([halfP, np.array([p])], axis=0)
    halfV = v


dt = 0.1
eulerP2 = np.array([np.array(P0)])
eulerV2 = np.array(V0)
for i in range(int(T / dt)):
    p = np.array(eulerP2[-1])
    v = np.array(eulerV2)
    a = F(p) / m
    v += (a * dt)
    p += v * dt

    if bounces and p[1] < 0:
        v *= coef_rest
        v[1] *= -1.0
        p[1] = max(0, p[1])
    eulerP2 = np.concatenate([eulerP2, np.array([p])], axis=0)
    eulerV2 = v

theoryP = np.array([P0])
for i in range(int(T / dt)):
    ti = i * dt
    theoryP = np.concatenate([theoryP, np.array([[P0[0] + V0[0] * ti, P0[1] + V0[1] * ti + 0.5 * -g * ti**2]])], axis=0)
    # theoryP = np.concatenate([theoryP, np.array([[P0[0] + V0[0] * ti + 0.5 * g[0] * ti**2, P0[1] + V0[1] * ti + 0.5 * g[1] * ti**2]])], axis=0)


# plt.plot(theoryP[:, 0], theoryP[:, 1], label = "Theory", linewidth = 2)
plt.plot(eulerP[:, 0], eulerP[:, 1], label = "Euler (dt = 1.0)")
plt.plot(semiP[:, 0], semiP[:, 1], label = "Semi-Euler (dt = 1.0)")
plt.plot(leapP[:, 0], leapP[:, 1], label = "Leapfrog (dt = 1.0)")
plt.plot(halfP[:, 0], halfP[:, 1], label = "Half (dt = 1.0)")
plt.plot(eulerP2[:, 0], eulerP2[:, 1], label = "Euler (dt = 0.1)")
plt.legend()
plt.show()
