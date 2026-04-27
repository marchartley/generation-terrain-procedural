import matplotlib.pyplot as plt
import numpy as np

m = 1.0
dt = 1.0
small_dt = 0.1
g = 1.0
T = 10.0
bounces = False
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


eulerP2 = np.array([np.array(P0)])
eulerV2 = np.array(V0)
for i in range(int(T / small_dt)):
    p = np.array(eulerP2[-1])
    v = np.array(eulerV2)
    a = F(p) / m
    v += (a * small_dt)
    p += v * small_dt

    if bounces and p[1] < 0:
        v *= coef_rest
        v[1] *= -1.0
        p[1] = max(0, p[1])
    eulerP2 = np.concatenate([eulerP2, np.array([p])], axis=0)
    eulerV2 = v

theoryP = np.array([P0])
for i in range(int(T / small_dt)):
    ti = i * small_dt
    theoryP = np.concatenate([theoryP, np.array([[P0[0] + V0[0] * ti, P0[1] + V0[1] * ti + 0.5 * -g * ti**2]])], axis=0)
    # theoryP = np.concatenate([theoryP, np.array([[P0[0] + V0[0] * ti + 0.5 * g[0] * ti**2, P0[1] + V0[1] * ti + 0.5 * g[1] * ti**2]])], axis=0)

fig, ax = plt.subplots()

lines = []
lines.append(ax.plot(theoryP[:, 0], theoryP[:, 1], label = "Theory", linewidth = 2)[0])
lines.append(ax.plot(eulerP[:, 0], eulerP[:, 1], label = f"Euler (dt = {dt})")[0])
lines.append(ax.plot(semiP[:, 0], semiP[:, 1], label = f"Semi-Euler (dt = {dt})")[0])
lines.append(ax.plot(leapP[:, 0], leapP[:, 1], label = f"Leapfrog (dt = {dt})")[0])
lines.append(ax.plot(halfP[:, 0], halfP[:, 1], label = f"Half (dt = {dt})")[0])
lines.append(ax.plot(eulerP2[:, 0], eulerP2[:, 1], label = f"Euler (dt = {small_dt})")[0])

leg = ax.legend()

lined = dict(zip(leg.get_lines() + leg.get_texts(), lines + lines))

for legline in leg.get_lines():
    legline.set_picker(True)
for text in leg.get_texts():
    text.set_picker(True)

def toggle(event):
    legline = event.artist
    line = lined[legline]
    line.set_visible(not line.get_visible())
    legline.set_alpha(1 if line.get_visible() else 0.2)
    ax.relim(visible_only=True)
    ax.autoscale_view()
    fig.canvas.draw_idle()

fig.canvas.mpl_connect('pick_event', toggle)

plt.show()
