import numpy as np


# Bezier Class representing a CUBIC bezier defined by four
# control points.
#
# at(t):            gets a point on the curve at t
# distance2(pt)      returns the closest distance^2 of
#                   pt and the curve
# closest(pt)       returns the point on the curve
#                   which is closest to pt
# maxes(pt)         plots the curve using matplotlib
class Bezier(object):
    exp3 = np.array([[3, 3], [2, 2], [1, 1], [0, 0]], dtype=np.float32)
    exp3_1 = np.array([[[3, 3], [2, 2], [1, 1], [0, 0]]], dtype=np.float32)
    exp4 = np.array([[4], [3], [2], [1], [0]], dtype=np.float32)
    boundaries = np.array([0, 1], dtype=np.float32)

    # Initialize the curve by assigning the control points.
    # Then create the coefficients.
    def __init__(self, points):
        assert isinstance(points, np.ndarray)
        assert points.dtype == np.float32
        self.points = points
        self.create_coefficients()

    # Create the coefficients of the bezier equation, bringing
    # the bezier in the form:
    # f(t) = a * t^3 + b * t^2 + c * t^1 + d
    #
    # The coefficients have the same dimensions as the control
    # points.
    def create_coefficients(self):
        points = self.points
        a = - points[0] + 3 * points[1] - 3 * points[2] + points[3]
        b = 3 * points[0] - 6 * points[1] + 3 * points[2]
        c = -3 * points[0] + 3 * points[1]
        d = points[0]
        self.coeffs = np.stack([a, b, c, d]).reshape(-1, 4, 2)

    # Return a point on the curve at the parameter t.
    def at(self, t):
        if type(t) != np.ndarray:
            t = np.array(t)
        pts = self.coeffs * np.power(t, self.exp3_1)
        return np.sum(pts, axis=1)

    # Return the closest DISTANCE (squared) between the point pt
    # and the curve.
    def distance2(self, pt):
        points, distances, index = self.measure_distance(pt)
        return distances[index]

    # Return the closest POINT between the point pt
    # and the curve.
    def closest(self, pt):
        points, distances, index = self.measure_distance(pt)
        return points[index]

    # Measure the distance^2 and closest point on the curve of
    # the point pt and the curve. This is done in a few steps:
    # 1     Define the distance^2 depending on the pt. I am
    #       using the squared distance because it is sufficient
    #       for comparing distances and doesn't have the over-
    #       head of an additional root operation.
    #       D(t) = (f(t) - pt)^2
    # 2     Get the roots of D'(t). These are the extremes of
    #       D(t) and contain the closest points on the unclipped
    #       curve. Only keep the minima by checking if
    #       D''(roots) > 0 and discard imaginary roots.
    # 3     Calculate the distances of the pt to the minima as
    #       well as the start and end of the curve and return
    #       the index of the shortest distance.
    #
    # This desmos graph is a helpful visualization.
    # https://www.desmos.com/calculator/ktglugn1ya
    def measure_distance(self, pt):
        coeffs = self.coeffs

        # These are the coefficients of the derivatives d/dx and d/(d/dx).
        da = 6 * np.sum(coeffs[0][0] * coeffs[0][0])
        db = 10 * np.sum(coeffs[0][0] * coeffs[0][1])
        dc = 4 * (np.sum(coeffs[0][1] * coeffs[0][1]) + 2 * np.sum(coeffs[0][0] * coeffs[0][2]))
        dd = 6 * (np.sum(coeffs[0][0] * (coeffs[0][3] - pt)) + np.sum(coeffs[0][1] * coeffs[0][2]))
        de = 2 * (np.sum(coeffs[0][2] * coeffs[0][2])) + 4 * np.sum(coeffs[0][1] * (coeffs[0][3] - pt))
        df = 2 * np.sum(coeffs[0][2] * (coeffs[0][3] - pt))

        dda = 5 * da
        ddb = 4 * db
        ddc = 3 * dc
        ddd = 2 * dd
        dde = de
        dcoeffs = np.stack([da, db, dc, dd, de, df])
        ddcoeffs = np.stack([dda, ddb, ddc, ddd, dde]).reshape(-1, 1)

        # Calculate the real extremes, by getting the roots of the first
        # derivativ of the distance function.
        extrema = np_real_roots(dcoeffs)
        # Remove the roots which are out of bounds of the clipped range [0, 1].
        # [future reference] https://stackoverflow.com/questions/47100903/deleting-every-3rd-element-of-a-tensor-in-tensorflow
        dd_clip = (np.sum(ddcoeffs * np.power(extrema, self.exp4)) >= 0) & (extrema > 0) & (extrema < 1)
        minima = extrema[dd_clip]

        # Add the start and end position as possible positions.
        potentials = np.concatenate((minima, self.boundaries))

        # Calculate the points at the possible parameters t and
        # get the index of the closest
        points = self.at(potentials.reshape(-1, 1, 1))
        distances = np.sum(np.square(points - pt), axis=1)
        index = np.argmin(distances)

        return points, distances, index

    # Point the curve to a matplotlib figure.
    # maxes         ... the axes of a matplotlib figure
    def plot(self, maxes):
        import matplotlib.path as mpath
        import matplotlib.patches as mpatches
        Path = mpath.Path
        pp1 = mpatches.PathPatch(
            Path(self.points, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
            fc="none")  # , transform=ax.transData)
        pp1.set_alpha(1)
        pp1.set_color('#00cc00')
        pp1.set_fill(False)
        pp2 = mpatches.PathPatch(
            Path(self.points, [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO]),
            fc="none")  # , transform=ax.transData)
        pp2.set_alpha(0.2)
        pp2.set_color('#666666')
        pp2.set_fill(False)

        maxes.scatter(*zip(*self.points), s=4, c=((0, 0.8, 1, 1), (0, 1, 0.5, 0.8), (0, 1, 0.5, 0.8),
                                                  (0, 0.8, 1, 1)))
        maxes.add_patch(pp2)
        maxes.add_patch(pp1)


# Wrapper around np.roots, but only returning real
# roots and ignoring imaginary results.
def np_real_roots(coefficients, EPSILON=1e-6):
    r = np.roots(coefficients)
    return r.real[abs(r.imag) < EPSILON]

import matplotlib.pyplot as plt
c = Bezier(np.array([[0, 0], [10, 5], [20, 5], [30, 30]], dtype=np.float32))

ax = plt.subplot()
c.plot(ax)
plt.show()
exit(0)








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
