import copy
import time
from dataclasses import dataclass
import os
from typing import List

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

# -----------------------------
# Utility: bilinear sampling
# -----------------------------
def _bilinear_sample(field, pts):
    """Bilinear sample a 2D field at continuous coordinates.
    field: HxW array
    pts: Nx2 array with (x,y) in pixel coords (x horizontal, y vertical)
    Returns: values (N,), and clamp mask (N,) saying which samples are fully inside.
    """
    H, W = field.shape
    x = pts[:, 0]
    y = pts[:, 1]

    x0 = np.floor(x).astype(int)
    y0 = np.floor(y).astype(int)
    x1 = x0 + 1
    y1 = y0 + 1

    # clamp for indexing; remember if inside
    inside = (x0 >= 0) & (x1 < W) & (y0 >= 0) & (y1 < H)

    x0c = np.clip(x0, 0, W-1); x1c = np.clip(x1, 0, W-1)
    y0c = np.clip(y0, 0, H-1); y1c = np.clip(y1, 0, H-1)

    Ia = field[y0c, x0c]
    Ib = field[y0c, x1c]
    Ic = field[y1c, x0c]
    Id = field[y1c, x1c]

    wa = (x1 - x) * (y1 - y)
    wb = (x - x0) * (y1 - y)
    wc = (x1 - x) * (y - y0)
    wd = (x - x0) * (y - y0)

    val = wa*Ia + wb*Ib + wc*Ic + wd*Id
    return val, inside

def _bilinear_sample_vec(field_x, field_y, pts):
    """Bilinear sample a 2D vector field at points pts (Nx2). Returns (vx, vy), inside mask."""
    vx, inside1 = _bilinear_sample(field_x, pts)
    vy, inside2 = _bilinear_sample(field_y, pts)
    return np.stack([vx, vy], axis=1), (inside1 & inside2)

# --------------------------------------
# Polyline helpers (closed/open curves)
# --------------------------------------
def roll(arr, shift):
    return np.roll(arr, shift=shift, axis=0)

def pairwise_diffs(x, closed):
    """Returns forward and backward neighbor differences for each vertex."""
    if closed:
        xm1 = roll(x, +1)
        xp1 = roll(x, -1)
    else:
        xm1 = np.vstack([x[0:1], x[:-1]])
        xp1 = np.vstack([x[1:], x[-1:]])
    return xp1 - x, x - xm1

def arc_lengths(x, closed):
    """Per-vertex arc-length weights (midpoint rule)."""
    if closed:
        xp1 = roll(x, -1)
        xm1 = roll(x, +1)
    else:
        xp1 = np.vstack([x[1:], x[-1:]])
        xm1 = np.vstack([x[0:1], x[:-1]])
    l = 0.5*(np.linalg.norm(xp1 - x, axis=1) + np.linalg.norm(x - xm1, axis=1))
    return l

def respaced_polyline(x, closed, n_points):
    """Resample polyline to equal arc-length spacing."""
    # compute cumulative arc length
    if closed:
        segs = roll(x, -1) - x
        Ls = np.linalg.norm(segs, axis=1)
        Lcum = np.concatenate([[0.0], np.cumsum(Ls)])
        Ltot = Lcum[-1]
        targets = np.linspace(0, Ltot, n_points+1)[:-1]  # exclude duplicate end
        # walk along
        new_pts = []
        j = 0
        for t in targets:
            while Lcum[j+1] < t:
                j += 1
            ratio = (t - Lcum[j]) / max(Ls[j], 1e-12)
            p = x[j] + ratio * (roll(x, -1)[j] - x[j])
            new_pts.append(p)
        return np.array(new_pts)
    else:
        segs = x[1:] - x[:-1]
        Ls = np.linalg.norm(segs, axis=1)
        Lcum = np.concatenate([[0.0], np.cumsum(Ls)])
        Ltot = Lcum[-1]
        targets = np.linspace(0, Ltot, n_points)
        new_pts = []
        j = 0
        for t in targets:
            while j < len(Ls)-1 and Lcum[j+1] < t:
                j += 1
            if Ls[j] < 1e-12:
                p = x[j].copy()
            else:
                ratio = (t - Lcum[j]) / Ls[j]
                p = x[j] + ratio * (x[j+1] - x[j])
            new_pts.append(p)
        return np.array(new_pts)

# --------------------------------------
# Main configuration dataclass
# --------------------------------------
@dataclass
class SnakeConfig:
    alpha: float = 0.05       # first-order regularization (tension)
    beta: float = 0.2         # second-order regularization (rigidity)
    gamma: float = 1.0        # external field weight (attraction to high f)
    eta: float = 0.5          # orientation (gradient alignment) weight
    lam_shape: float = 0.0    # length target weight
    L_target: float = None    # target length in pixels (if None, disabled)
    tau: float = 0.1          # step size
    max_iter: int = 400
    closed: bool = True
    freeze_grad: bool = True  # treat grad f as constant per iteration
    normalize_grad: bool = True  # use gradient direction only for orientation
    respacing_every: int = 20
    n_points: int = 200       # number of vertices
    clamp_to_image: bool = True  # keep vertices inside image bounds
    report_every: int = 50
    fixed_vertices: List[int] = None

    # --- Chan–Vese region term parameters ---
    use_chanvese: bool = False
    cv_lambda1: float = 1.0   # weight for inside term
    cv_lambda2: float = 1.0   # weight for outside term
    cv_nu: float = 0.0        # pure area bias term (ν)
    cv_weight: float = 1.0    # global scaling for CV term

    # --- Area constraint (target polygon area) ---
    lam_area: float = 0.0     # weight of (|A| - A_target)^2
    A_target: float = None    # target area in pixels^2 (if None, disabled)

# --------------------------------------
# Snake optimizer
# --------------------------------------
class Snake:
    def __init__(self, field, grad=None, cfg: SnakeConfig = SnakeConfig()):
        """
        field: HxW float32 in [0,1], larger = more attractive (E_ext = -gamma * ∫ f ds)
        grad: optional precomputed gradient as (gy, gx) arrays same shape as field
        """
        self.f = field.astype(np.float32)
        self.H, self.W = self.f.shape
        if grad is None:
            gy, gx = np.gradient(self.f)
            self.gx = gx.astype(np.float32)
            self.gy = gy.astype(np.float32)
        else:
            self.gy, self.gx = grad
        self.cfg = cfg

        # Precompute grid of pixel centers for Chan–Vese region statistics
        yy, xx = np.mgrid[0:self.H, 0:self.W]
        self._cv_grid_pts = np.stack([xx.ravel() + 0.5, yy.ravel() + 0.5], axis=1)

    def _sample_f_and_grad(self, pts):
        vals, inside = _bilinear_sample(self.f, pts)
        g, inside2 = _bilinear_sample_vec(self.gx, self.gy, pts)  # note order (x,y)
        inside = inside & inside2
        return vals, g, inside

    def _internal_gradient(self, x):
        """Classic snake internal forces: 1st and 2nd finite differences (periodic if closed)."""
        if self.cfg.closed:
            xm1 = roll(x, +1); xp1 = roll(x, -1)
            xm2 = roll(x, +2); xp2 = roll(x, -2)
        else:
            xm1 = np.vstack([x[0:1], x[:-1]])
            xp1 = np.vstack([x[1:], x[-1:]])
            xm2 = np.vstack([x[0:1], x[:-2], x[-2:-1]])
            xp2 = np.vstack([x[2:], x[-1:], x[-1:]])
        alpha = self.cfg.alpha
        beta = self.cfg.beta
        g1 = 2*alpha*(2*x - xm1 - xp1)
        g2 = 2*beta*(6*x - 4*(xm1 + xp1) + (xm2 + xp2))
        return g1 + g2

    def _external_gradient(self, x):
        # E_ext = -gamma * sum f(x_i) * l_i  -> grad = -gamma * l_i * grad f(x_i)
        vals, g, inside = self._sample_f_and_grad(x)
        l = arc_lengths(x, self.cfg.closed)[:, None]
        return -self.cfg.gamma * l * g

    def _length_target_gradient(self, x):
        if self.cfg.lam_shape <= 0.0 or self.cfg.L_target is None:
            return np.zeros_like(x)
        if self.cfg.closed:
            xp1 = roll(x, -1); xm1 = roll(x, +1)
        else:
            xp1 = np.vstack([x[1:], x[-1:]])
            xm1 = np.vstack([x[0:1], x[:-1]])
        seg_fwd = xp1 - x
        seg_bwd = x - xm1
        L = np.sum(np.linalg.norm(seg_fwd, axis=1))
        e1 = seg_bwd / (np.linalg.norm(seg_bwd, axis=1, keepdims=True) + 1e-8)
        e2 = seg_fwd / (np.linalg.norm(seg_fwd, axis=1, keepdims=True) + 1e-8)
        dLdx = e1 - e2
        return -2.0 * self.cfg.lam_shape * (self.cfg.L_target - L) * dLdx

    def _orientation_force(self, x):
        """Heuristic 'frozen-gradient alignment' force that encourages tangent || grad f."""
        # Tangent using central difference
        if self.cfg.closed:
            xp1 = roll(x, -1); xm1 = roll(x, +1)
        else:
            xp1 = np.vstack([x[1:], x[-1:]])
            xm1 = np.vstack([x[0:1], x[:-1]])
        t = xp1 - xm1
        t_norm = np.linalg.norm(t, axis=1, keepdims=True) + 1e-8
        t_hat = t / t_norm

        # Sample gradient
        _, g, _ = self._sample_f_and_grad(x)
        if self.cfg.normalize_grad:
            g_norm = np.linalg.norm(g, axis=1, keepdims=True) + 1e-8
            g_hat = g / g_norm
        else:
            g_hat = g

        # Misalignment component of tangent orthogonal to gradient direction
        dot = np.sum(t_hat * g_hat, axis=1, keepdims=True)
        m = t_hat - np.maximum(dot, 0.0) * g_hat

        # Convert this "directional mismatch" into vertex forces using a divergence-like stencil
        if self.cfg.closed:
            mp1 = roll(m, -1)
            mm1 = roll(m, +1)
        else:
            mp1 = np.vstack([m[1:], m[-1:]])
            mm1 = np.vstack([m[0:1], m[:-1]])

        F = self.cfg.eta * (mp1 - mm1)

        # Weight by local arc-length to be invariant to sampling density
        l = arc_lengths(x, self.cfg.closed)[:, None]
        return F * l

    # --- NEW: signed polygon area (positive if vertices are CCW) ---
    def _polygon_area(self, x):
        """Signed area of closed polygon x (Nx2)."""
        xp1 = roll(x, -1)
        return 0.5 * np.sum(x[:, 0] * xp1[:, 1] - xp1[:, 0] * x[:, 1])

    # --- NEW: area constraint gradient ---
    def _area_constraint_gradient(self, x):
        """
        E_area = lam_area * (|A| - A_target)^2
        ∂A/∂x_i = 0.5 * [ y_{i+1} - y_{i-1}, x_{i-1} - x_{i+1} ]
        """
        cfg = self.cfg
        if (not cfg.closed) or cfg.lam_area <= 0.0 or cfg.A_target is None:
            return np.zeros_like(x)

        xp1 = roll(x, -1)
        xm1 = roll(x, +1)

        # gradient of signed area
        dA_dx = 0.5 * np.stack([xp1[:, 1] - xm1[:, 1],
                                xm1[:, 0] - xp1[:, 0]], axis=1)

        A = self._polygon_area(x)
        A_abs = abs(A)
        # derivative of |A| is sign(A)*dA_dx (for A != 0)
        signA = 1.0 if A >= 0.0 else -1.0
        dAbsA_dx = signA * dA_dx

        # E = lam * (|A| - A_target)^2
        # ∂E/∂x = 2 * lam * (|A| - A_target) * ∂|A|/∂x
        return 2.0 * cfg.lam_area * (A_abs - cfg.A_target) * dAbsA_dx

    # --- NEW: Chan–Vese region-based term ---
    def _chanvese_region_gradient(self, x):
        """
        Region-based force à la Chan–Vese:

        E_CV ∝ λ1 ∫_inside (I - c1)^2 + λ2 ∫_outside (I - c2)^2 + ν * Area(inside)

        We approximate gradient descent by pushing vertices along the normal with
        speed:
            F = ν + λ1 (I(x) - c1)^2 - λ2 (I(x) - c2)^2
        and then:
            grad_E ≈ F * n_hat * (arc-length weight)
        """
        cfg = self.cfg
        if (not cfg.closed) or (not cfg.use_chanvese) or cfg.cv_weight == 0.0:
            return np.zeros_like(x)

        from matplotlib.path import Path

        H, W = self.H, self.W
        I = self.f

        # 1) build polygon path and compute inside/outside mask
        path = Path(x)  # closed polygon
        inside_flat = path.contains_points(self._cv_grid_pts)
        inside = inside_flat.reshape(H, W)

        n_in = inside.sum()
        n_total = H * W
        n_out = n_total - n_in

        # guard against degenerate cases
        if n_in == 0 or n_out == 0:
            return np.zeros_like(x)

        # 2) region means c1 (inside), c2 (outside)
        c1 = float(I[inside].mean())
        c2 = float(I[~inside].mean())

        # 3) intensity sampled along the curve
        I_curve, _, _ = self._sample_f_and_grad(x)  # values at vertices, shape (N,)

        lam1 = cfg.cv_lambda1
        lam2 = cfg.cv_lambda2
        nu   = cfg.cv_nu

        # Chan–Vese “speed” term along the normal (scalar per vertex)
        F = nu + lam1 * (I_curve - c1) ** 2 - lam2 * (I_curve - c2) ** 2  # (N,)

        # 4) compute unit normals
        xp1 = roll(x, -1)
        xm1 = roll(x, +1)
        t = xp1 - xm1  # approximate tangent
        # rotate tangent to get normal: n = (t_y, -t_x)
        n = np.stack([t[:, 1], -t[:, 0]], axis=1)
        n_norm = np.linalg.norm(n, axis=1, keepdims=True) + 1e-8
        n_hat = n / n_norm

        # 5) weight by local arc-length to make it invariant to sampling
        l = arc_lengths(x, cfg.closed)[:, None]

        # gradient of energy (same shape as x)
        g_cv = cfg.cv_weight * F[:, None] * n_hat * l
        return g_cv.astype(np.float32)

    def _clamp(self, x):
        if not self.cfg.clamp_to_image:
            return x
        x[:, 0] = np.clip(x[:, 0], 0, self.W - 1)
        x[:, 1] = np.clip(x[:, 1], 0, self.H - 1)
        return x

    def optimize(self, x0, verbose=False):
        paths = [x0.copy()]
        x = x0.copy().astype(np.float32)
        cfg = self.cfg
        for it in range(cfg.max_iter):
            g_int = self._internal_gradient(x)
            g_ext = self._external_gradient(x)
            g_len = self._length_target_gradient(x)
            g_ori = -self._orientation_force(x)

            # Chan–Vese region term and area constraint
            g_cv   = self._chanvese_region_gradient(x)
            g_area = self._area_constraint_gradient(x)

            # total gradient
            g = g_int + g_ext + g_len + g_ori + g_cv + g_area

            x = x - cfg.tau * g
            x = self._clamp(x)

            if self.cfg.fixed_vertices is not None:
                x[self.cfg.fixed_vertices] = x0[self.cfg.fixed_vertices]

            if cfg.respacing_every > 0 and (it+1) % cfg.respacing_every == 0:
                x = respaced_polyline(x, cfg.closed, cfg.n_points)
            paths.append(x.copy())

            if verbose and (it % max(1,cfg.report_every) == 0):
                # crude energy report
                E_int = 0.5*np.sum((x - roll(x,+1))**2) + 0.5*np.sum((x - 2*roll(x,+1) + roll(x,+2))**2)
                vals, _, _ = self._sample_f_and_grad(x)
                l = arc_lengths(x, cfg.closed)
                E_ext = -cfg.gamma * np.sum(vals * l)
                print(f"iter {it:4d}: E_int~{E_int:.2f}  E_ext~{E_ext:.2f}")
        return x, paths

# --------------------------------------
# Initialization helpers
# --------------------------------------
def init_circle(W, H, n_points=200, margin_frac=0.35, center=None):
    if center is None:
        cx, cy = W/2.0, H/2.0
    else:
        cx, cy = center
    r = (1 - margin_frac) * min(W, H) / 2
    t = np.linspace(0, 2*np.pi, n_points, endpoint=False)
    x = np.stack([cx + r*np.cos(t), cy + r*np.sin(t)], axis=1)
    return x.astype(np.float32)

def init_rectangle(W, H, n_points=200, margin_frac=0.2):
    # Axis-aligned rectangle approximated by equal steps
    m = margin_frac
    pts = np.array([
        [m*W, m*H],
        [W-m*W, m*H],
        [W-m*W, H-m*H],
        [m*W, H-m*H],
    ], dtype=np.float32)
    # resample
    return respaced_polyline(np.vstack([pts, pts[0]]), True, n_points)


# -----------------------------
# Load image & set up snake
# -----------------------------
img_path = "/home/marc/generation-terrain-procedural/EnvObjRendering/2025-12-04__14-28-37/material_polyp.png"

img = Image.open(img_path).convert("I;16")
img = img.resize((50, 50))
f = np.asarray(img, dtype=np.float32) / 65535.0

H, W = f.shape

# --- Config ---
cfg = SnakeConfig(
    alpha=0.03,         # first-order regularization (tension)
    beta=0.05,         # second-order regularization (rigidity)
    gamma=-2.0,        # external field weight (attraction to high f)
    eta=0.0,          # orientation (gradient alignment) weight
    lam_shape=0.03,     # length target weight
    L_target=W,        # target length in pixels (if None, disabled)
    tau=0.1,          # step size
    max_iter=1000,
    closed=True,
    freeze_grad=True,
    normalize_grad=False,
    respacing_every=15,
    n_points=30,
    clamp_to_image=True,
    report_every=100,
    fixed_vertices=[0],

    # --- Chan–Vese usage (tune these) ---
    use_chanvese=True,   # set True to enable CV term
    cv_lambda1=1.0,
    cv_lambda2=1.0,
    cv_nu=0.0,
    cv_weight=5.0,

    # --- Area constraint (optional) ---
    lam_area=0.01,         # set >0 to constrain area
    A_target=200.0
)

# --- Initialize base curve ---
base_x0 = init_rectangle(W, H, n_points=cfg.n_points, margin_frac=.1)
# base_x0 = init_circle(W, H, n_points=cfg.n_points, margin_frac=.1)

# If you want area constraint, initialize A_target to the initial area:
# snk_tmp = Snake(f, cfg=cfg)
# cfg.A_target = abs(snk_tmp._polygon_area(base_x0))
# cfg.lam_area = 1e-4  # for example

# --- Optimizer (reused across updates) ---
snk = Snake(f, cfg=cfg)

# --- Matplotlib setup ---
fig, ax = plt.subplots(figsize=(7, 5))
gy, gx = np.gradient(f.astype(float))
grad_mag = np.sqrt(gx**2 + gy**2)
H, W = f.shape
Y, X = np.mgrid[0:H, 0:W]

# draw isolines
ax.contour(X, Y, f, levels=20, cmap="gray")  # 20 contour levels
ax.axis("off")

from matplotlib.collections import LineCollection
paths_coll = LineCollection([], linewidths=1, alpha=0.6)  # thin, semi-transparent
ax.add_collection(paths_coll)

import matplotlib as mpl
cmap = plt.get_cmap("viridis")
(line,) = ax.plot([], [], lw=4.5)
(draw_line,) = ax.plot([], [], linestyle="--", linewidth=1)

DEBOUNCE = 0.12
_last_run = 0.0
_last_xy = (None, None)
_lastPath = None

#  state for drawing
is_drawing = False
draw_pts = []   # list of (x,y) while dragging

def run_snake_at(x, y, currentPath=None, useMouse=True):
    """Run the snake with the first point forced at (x, y) using current base_x0."""
    if currentPath is None:
        x0 = base_x0.copy()
    else:
        x0 = currentPath.copy()
    if useMouse:
        x0[0] = [x, y]
    # Optimize (your optimize returns (result_path, paths))
    prevFreeze = snk.cfg.fixed_vertices
    if not useMouse:
        snk.cfg.fixed_vertices = []
    result_path, paths = snk.optimize(x0, verbose=False)
    snk.cfg.fixed_vertices = prevFreeze
    xT = paths[-1]
    return xT, paths

def on_key_press(event):
    if event.key == 'enter':
        global _last_run, _last_xy, draw_pts, line, draw_line
        xT, allPaths = run_snake_at(*_last_xy, currentPath=None, useMouse=False)
        # Build segments for the LineCollection (one polyline per iteration)
        if cfg.closed:
            segs = [np.vstack([p, p[0]]) for p in allPaths[50::10]]   # close each loop
        else:
            segs = [p for p in allPaths[50::10]]                      # open polylines

        line.set_data([], [])
        draw_line.set_data([], [])
        n = len(segs)
        for iSnake in range(n):
            paths_coll.set_segments([segs[iSnake]])
            colors = "blue"
            paths_coll.set_color(colors)
            paths_coll.set_alpha(1.0)
            paths_coll.set_linewidth(4.5)

            fig.canvas.draw_idle()
            fig.savefig(f"snake-test-{iSnake}.png")
            print(f"Saved {iSnake}/{n}")
        paths_coll.set_linewidth(1.0)
        paths_coll.set_alpha(0.6)

def on_mouse_move(event):
    global _last_run, _last_xy, draw_pts
    if event.inaxes != ax or event.xdata is None or event.ydata is None:
        return
    x, y = float(event.xdata), float(event.ydata)

    if is_drawing:
        # add points as the mouse moves to build the initial curve
        if not draw_pts or (abs(x - draw_pts[-1][0]) + abs(y - draw_pts[-1][1])) > 0.75:
            draw_pts.append((x, y))
            xs, ys = zip(*draw_pts)
            draw_line.set_data(xs, ys)
            fig.canvas.draw_idle()
        return

    # normal hover: throttle recomputes
    t = time.monotonic()
    if (t - _last_run) < DEBOUNCE:
        return

    lx, ly = _last_xy
    if lx is not None and (abs(x - lx) < 1.0 and abs(y - ly) < 1.0):
        return

    _last_run = t
    _last_xy = (x, y)

    xT, allPaths = run_snake_at(x, y, currentPath=None)

    # Build segments for the LineCollection (one polyline per iteration)
    if cfg.closed:
        segs = [np.vstack([p, p[0]]) for p in allPaths[50::50]]   # close each loop
    else:
        segs = [p for p in allPaths[50::50]]                      # open polylines

    # Update the collection with all paths
    paths_coll.set_segments(segs)

    # Optional: color by iteration (earlier = lighter, later = darker)
    n = len(segs)
    colors = [cmap(i/(n-1 if n>1 else 1)) for i in range(n)]
    paths_coll.set_color(colors)

    # Keep the final path emphasized with the existing Line2D
    if cfg.closed:
        line.set_data(np.r_[xT[:, 0], xT[0, 0]], np.r_[xT[:, 1], xT[0, 1]])
    else:
        line.set_data(xT[:, 0], xT[:, 1])

    fig.canvas.draw_idle()

def on_button_press(event):
    # begin drawing a custom initial curve
    global is_drawing, draw_pts
    if event.inaxes != ax or event.xdata is None or event.ydata is None:
        return
    is_drawing = True    # start drawing
    draw_pts = [(float(event.xdata), float(event.ydata))]
    draw_line.set_data([event.xdata], [event.ydata])
    fig.canvas.draw_idle()

def on_button_release(event):
    # finalize the drawn curve, resample to cfg.n_points, set as new base_x0
    global is_drawing, base_x0, draw_pts
    if not is_drawing:
        return
    is_drawing = False

    if len(draw_pts) >= 2:
        pts = np.array(draw_pts, dtype=np.float32)
        # If closed, close the loop if the endpoints are near
        if cfg.closed:
            if np.linalg.norm(pts[0] - pts[-1]) > 1.0:
                pts = np.vstack([pts, pts[0]])
        # resample to uniform spacing with the same open/closed setting
        base_x0 = respaced_polyline(pts, cfg.closed, cfg.n_points)
        cfg.L_target = np.sum(arc_lengths(base_x0, cfg.closed)) * 1.1
        # update preview line to show the new base curve
        if cfg.closed:
            draw_line.set_data(np.r_[base_x0[:,0], base_x0[0,0]], np.r_[base_x0[:,1], base_x0[0,1]])
        else:
            draw_line.set_data(base_x0[:,0], base_x0[:,1])
    else:
        # not enough points: keep previous base_x0
        draw_line.set_data([], [])

    draw_pts = []
    fig.canvas.draw_idle()

# connect events
cid_move   = fig.canvas.mpl_connect('motion_notify_event', on_mouse_move)
cid_press  = fig.canvas.mpl_connect('button_press_event', on_button_press)
cid_release= fig.canvas.mpl_connect('button_release_event', on_button_release)
cid_key    = fig.canvas.mpl_connect('key_press_event', on_key_press)

plt.show()
exit()















import copy
import time
from dataclasses import dataclass
import os
from typing import List

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

# -----------------------------
# Utility: bilinear sampling
# -----------------------------
def _bilinear_sample(field, pts):
    """Bilinear sample a 2D field at continuous coordinates.
    field: HxW array
    pts: Nx2 array with (x,y) in pixel coords (x horizontal, y vertical)
    Returns: values (N,), and clamp mask (N,) saying which samples are fully inside.
    """
    H, W = field.shape
    x = pts[:, 0]
    y = pts[:, 1]

    x0 = np.floor(x).astype(int)
    y0 = np.floor(y).astype(int)
    x1 = x0 + 1
    y1 = y0 + 1

    # clamp for indexing; remember if inside
    inside = (x0 >= 0) & (x1 < W) & (y0 >= 0) & (y1 < H)

    x0c = np.clip(x0, 0, W-1); x1c = np.clip(x1, 0, W-1)
    y0c = np.clip(y0, 0, H-1); y1c = np.clip(y1, 0, H-1)

    Ia = field[y0c, x0c]
    Ib = field[y0c, x1c]
    Ic = field[y1c, x0c]
    Id = field[y1c, x1c]

    wa = (x1 - x) * (y1 - y)
    wb = (x - x0) * (y1 - y)
    wc = (x1 - x) * (y - y0)
    wd = (x - x0) * (y - y0)

    val = wa*Ia + wb*Ib + wc*Ic + wd*Id
    return val, inside

def _bilinear_sample_vec(field_x, field_y, pts):
    """Bilinear sample a 2D vector field at points pts (Nx2). Returns (vx, vy), inside mask."""
    vx, inside1 = _bilinear_sample(field_x, pts)
    vy, inside2 = _bilinear_sample(field_y, pts)
    return np.stack([vx, vy], axis=1), (inside1 & inside2)

# --------------------------------------
# Polyline helpers (closed/open curves)
# --------------------------------------
def roll(arr, shift):
    return np.roll(arr, shift=shift, axis=0)

def pairwise_diffs(x, closed):
    """Returns forward and backward neighbor differences for each vertex."""
    if closed:
        xm1 = roll(x, +1)
        xp1 = roll(x, -1)
    else:
        xm1 = np.vstack([x[0:1], x[:-1]])
        xp1 = np.vstack([x[1:], x[-1:]])
    return xp1 - x, x - xm1

def arc_lengths(x, closed):
    """Per-vertex arc-length weights (midpoint rule)."""
    if closed:
        xp1 = roll(x, -1)
        xm1 = roll(x, +1)
    else:
        xp1 = np.vstack([x[1:], x[-1:]])
        xm1 = np.vstack([x[0:1], x[:-1]])
    l = 0.5*(np.linalg.norm(xp1 - x, axis=1) + np.linalg.norm(x - xm1, axis=1))
    return l

def respaced_polyline(x, closed, n_points):
    """Resample polyline to equal arc-length spacing."""
    # compute cumulative arc length
    if closed:
        segs = roll(x, -1) - x
        Ls = np.linalg.norm(segs, axis=1)
        Lcum = np.concatenate([[0.0], np.cumsum(Ls)])
        Ltot = Lcum[-1]
        targets = np.linspace(0, Ltot, n_points+1)[:-1]  # exclude duplicate end
        # walk along
        new_pts = []
        j = 0
        for t in targets:
            while Lcum[j+1] < t:
                j += 1
            ratio = (t - Lcum[j]) / max(Ls[j], 1e-12)
            p = x[j] + ratio * (roll(x, -1)[j] - x[j])
            new_pts.append(p)
        return np.array(new_pts)
    else:
        segs = x[1:] - x[:-1]
        Ls = np.linalg.norm(segs, axis=1)
        Lcum = np.concatenate([[0.0], np.cumsum(Ls)])
        Ltot = Lcum[-1]
        targets = np.linspace(0, Ltot, n_points)
        new_pts = []
        j = 0
        for t in targets:
            while j < len(Ls)-1 and Lcum[j+1] < t:
                j += 1
            if Ls[j] < 1e-12:
                p = x[j].copy()
            else:
                ratio = (t - Lcum[j]) / Ls[j]
                p = x[j] + ratio * (x[j+1] - x[j])
            new_pts.append(p)
        return np.array(new_pts)

# --------------------------------------
# Main configuration dataclass
# --------------------------------------
@dataclass
class SnakeConfig:
    alpha: float = 0.05       # first-order regularization (tension)
    beta: float = 0.2         # second-order regularization (rigidity)
    gamma: float = 1.0        # external field weight (attraction to high f)
    eta: float = 0.5          # orientation (gradient alignment) weight
    lam_shape: float = 0.0    # length target weight
    L_target: float = None    # target length in pixels (if None, disabled)
    tau: float = 0.1          # step size
    max_iter: int = 400
    closed: bool = True
    freeze_grad: bool = True  # treat grad f as constant per iteration
    normalize_grad: bool = True  # use gradient direction only for orientation
    respacing_every: int = 20
    n_points: int = 200       # number of vertices
    clamp_to_image: bool = True  # keep vertices inside image bounds
    report_every: int = 50
    fixed_vertices: List[int] = None

# --------------------------------------
# Snake optimizer
# --------------------------------------
class Snake:
    def __init__(self, field, grad=None, cfg: SnakeConfig = SnakeConfig()):
        """
        field: HxW float32 in [0,1], larger = more attractive (E_ext = -gamma * ∫ f ds)
        grad: optional precomputed gradient as (gy, gx) arrays same shape as field
        """
        self.f = field.astype(np.float32)
        self.H, self.W = self.f.shape
        if grad is None:
            gy, gx = np.gradient(self.f)
            self.gx = gx.astype(np.float32)
            self.gy = gy.astype(np.float32)
        else:
            self.gy, self.gx = grad
        self.cfg = cfg

    def _sample_f_and_grad(self, pts):
        vals, inside = _bilinear_sample(self.f, pts)
        g, inside2 = _bilinear_sample_vec(self.gx, self.gy, pts)  # note order (x,y)
        inside = inside & inside2
        return vals, g, inside

    def _internal_gradient(self, x):
        """Classic snake internal forces: 1st and 2nd finite differences (periodic if closed)."""
        if self.cfg.closed:
            xm1 = roll(x, +1); xp1 = roll(x, -1)
            xm2 = roll(x, +2); xp2 = roll(x, -2)
        else:
            xm1 = np.vstack([x[0:1], x[:-1]])
            xp1 = np.vstack([x[1:], x[-1:]])
            xm2 = np.vstack([x[0:1], x[:-2], x[-2:-1]])
            xp2 = np.vstack([x[2:], x[-1:], x[-1:]])
        alpha = self.cfg.alpha
        beta = self.cfg.beta
        g1 = 2*alpha*(2*x - xm1 - xp1)
        g2 = 2*beta*(6*x - 4*(xm1 + xp1) + (xm2 + xp2))
        return g1 + g2

    def _external_gradient(self, x):
        # E_ext = -gamma * sum f(x_i) * l_i  -> grad = -gamma * l_i * grad f(x_i)
        vals, g, inside = self._sample_f_and_grad(x)
        l = arc_lengths(x, self.cfg.closed)[:, None]
        return -self.cfg.gamma * l * g

    def _length_target_gradient(self, x):
        if self.cfg.lam_shape <= 0.0 or self.cfg.L_target is None:
            return np.zeros_like(x)
        if self.cfg.closed:
            xp1 = roll(x, -1); xm1 = roll(x, +1)
        else:
            xp1 = np.vstack([x[1:], x[-1:]])
            xm1 = np.vstack([x[0:1], x[:-1]])
        seg_fwd = xp1 - x
        seg_bwd = x - xm1
        L = np.sum(np.linalg.norm(seg_fwd, axis=1))
        e1 = seg_bwd / (np.linalg.norm(seg_bwd, axis=1, keepdims=True) + 1e-8)
        e2 = seg_fwd / (np.linalg.norm(seg_fwd, axis=1, keepdims=True) + 1e-8)
        dLdx = e1 - e2
        return -2.0 * self.cfg.lam_shape * (self.cfg.L_target - L) * dLdx

    def _orientation_force(self, x):
        """Heuristic 'frozen-gradient alignment' force that encourages tangent || grad f.
        We minimize ||t_hat x g_hat||^2 by reducing the component of the tangent
        orthogonal to the gradient direction (no Hessians; gradient treated as constant).
        Returns a force with the same units as a gradient of energy.
        """
        # if self.cfg.eta <= 0.0:
        #     return np.zeros_like(x)

        # Tangent using central difference
        if self.cfg.closed:
            xp1 = roll(x, -1); xm1 = roll(x, +1)
        else:
            xp1 = np.vstack([x[1:], x[-1:]])
            xm1 = np.vstack([x[0:1], x[:-1]])
        t = xp1 - xm1
        t_norm = np.linalg.norm(t, axis=1, keepdims=True) + 1e-8
        t_hat = t / t_norm

        # Sample gradient
        _, g, _ = self._sample_f_and_grad(x)
        if self.cfg.normalize_grad:
            g_norm = np.linalg.norm(g, axis=1, keepdims=True) + 1e-8
            g_hat = g / g_norm
        else:
            g_hat = g

        # Misalignment component of tangent orthogonal to gradient direction
        # m = t_hat - (t_hat·g_hat) g_hat
        dot = np.sum(t_hat * g_hat, axis=1, keepdims=True)
        # dot = np.clip(dot, 0, 1)
        # m = t_hat - dot * g_hat
        m = t_hat - np.maximum(dot, 0.0) * g_hat

        # Convert this "directional mismatch" into vertex forces using a divergence-like stencil
        # so that reducing m smooths along the curve while aligning it to g_hat.
        if self.cfg.closed:
            mp1 = roll(m, -1)
            mm1 = roll(m, +1)
        else:
            mp1 = np.vstack([m[1:], m[-1:]])
            mm1 = np.vstack([m[0:1], m[:-1]])

            # mp1 = np.clip(mp1, 0, 1)
            # mm1 = np.clip(mm1, 0, 1)

        # Force tries to reduce local change in m along the curve:
        F = self.cfg.eta * (mp1 - mm1)

        # Weight by local arc-length to be invariant to sampling density
        l = arc_lengths(x, self.cfg.closed)[:, None]
        return F * l

    def _clamp(self, x):
        if not self.cfg.clamp_to_image:
            return x
        x[:, 0] = np.clip(x[:, 0], 0, self.W - 1)
        x[:, 1] = np.clip(x[:, 1], 0, self.H - 1)
        return x

    def optimize(self, x0, verbose=False):
        paths = [x0.copy()]
        x = x0.copy().astype(np.float32)
        cfg = self.cfg
        for it in range(cfg.max_iter):
            g_int = self._internal_gradient(x)
            g_ext = self._external_gradient(x)
            g_len = self._length_target_gradient(x)
            g_ori = -self._orientation_force(x)

            g = g_int + g_ext + g_len + g_ori
            x = x - cfg.tau * g
            x = self._clamp(x)

            if self.cfg.fixed_vertices is not None:
                x[self.cfg.fixed_vertices] = x0[self.cfg.fixed_vertices]

            if cfg.respacing_every > 0 and (it+1) % cfg.respacing_every == 0:
                x = respaced_polyline(x, cfg.closed, cfg.n_points)
            paths.append(x.copy())

            if verbose and (it % max(1,cfg.report_every) == 0):
                # crude energy report
                E_int = 0.5*np.sum((x - roll(x,+1))**2) + 0.5*np.sum((x - 2*roll(x,+1) + roll(x,+2))**2)
                vals, _, _ = self._sample_f_and_grad(x)
                l = arc_lengths(x, cfg.closed)
                E_ext = -cfg.gamma * np.sum(vals * l)
                print(f"iter {it:4d}: E_int~{E_int:.2f}  E_ext~{E_ext:.2f}")
        return x, paths

# --------------------------------------
# Initialization helpers
# --------------------------------------
def init_circle(W, H, n_points=200, margin_frac=0.35, center=None):
    if center is None:
        cx, cy = W/2.0, H/2.0
    else:
        cx, cy = center
    r = (1 - margin_frac) * min(W, H) / 2
    t = np.linspace(0, 2*np.pi, n_points, endpoint=False)
    x = np.stack([cx + r*np.cos(t), cy + r*np.sin(t)], axis=1)
    return x.astype(np.float32)

def init_rectangle(W, H, n_points=200, margin_frac=0.2):
    # Axis-aligned rectangle approximated by equal steps
    m = margin_frac
    pts = np.array([
        [m*W, m*H],
        [W-m*W, m*H],
        [W-m*W, H-m*H],
        [m*W, H-m*H],
    ], dtype=np.float32)
    # resample
    return respaced_polyline(np.vstack([pts, pts[0]]), True, n_points)



# img_path = "/home/marc/generation-terrain-procedural/EnvObjRendering/2025-12-04__13-24-43/material_polyp.png" #"snake_test.png"
img_path = "/home/marc/generation-terrain-procedural/EnvObjRendering/2025-12-04__14-28-37/material_polyp.png" #"snake_test.png"
# img = Image.open(img_path).convert("L")
# img = img.resize((512, 512))
# f = np.asarray(img, dtype=np.float32) / 255.0
img = Image.open(img_path).convert("I;16")
f = np.asarray(img, dtype=np.float32) / 65535.0
#
# plt.imshow(img)
# plt.show()

H, W = f.shape

# --- Config ---
cfg = SnakeConfig(
    alpha=0.3,         # first-order regularization (tension)
    beta=0.05,          # second-order regularization (rigidity)
    gamma=10.0,          # external field weight (attraction to high f)
    eta=-1.0,           # orientation (gradient alignment) weight
    lam_shape=0.3,    # length target weight
    L_target= W,        # target length in pixels (if None, disabled)
    tau=0.01,           # step size
    max_iter=1000,
    closed=True,
    freeze_grad=True,
    normalize_grad=True,
    respacing_every=15,
    n_points=100,
    clamp_to_image=True,
    report_every=100,
    fixed_vertices=[0]
)

# --- Initialize base curve (will reinitialize on each mouse move) ---
# base_x0 = init_circle(W, H, n_points=cfg.n_points, radius_frac=0.35)
# base_x0 = init_rectangle(W, H, n_points=cfg.n_points, margin_frac=.1)
base_x0 = init_circle(W, H, n_points=cfg.n_points, margin_frac=.1)

# --- Optimizer (reused across updates) ---
snk = Snake(f, cfg=cfg)

# --- Matplotlib setup ---
fig, ax = plt.subplots(figsize=(7, 5))
gy, gx = np.gradient(f.astype(float))
grad_mag = np.sqrt(gx**2 + gy**2)
H, W = f.shape
Y, X = np.mgrid[0:H, 0:W]

# draw isolines
ax.contour(X, Y, f, levels=20, cmap="gray")  # 20 contour levels
# ax.imshow(f, cmap="gray", origin="upper")
# ax.imshow(grad_mag, cmap="gray", origin="upper")
ax.axis("off")

# single Line2D we’ll update
from matplotlib.collections import LineCollection

# after you create fig, ax ...
paths_coll = LineCollection([], linewidths=1, alpha=0.6)  # thin, semi-transparent
ax.add_collection(paths_coll)

# optional: color progression from early→late
import matplotlib as mpl
cmap = plt.get_cmap("viridis")
(line,) = ax.plot([], [], lw=4.5)
(draw_line,) = ax.plot([], [], linestyle="--", linewidth=1)

# simple throttling so we don’t fire too often (in seconds)
DEBOUNCE = 0.12
_last_run = 0.0
_last_xy = (None, None)
_lastPath = None

#  state for drawing
is_drawing = False
draw_pts = []   # list of (x,y) while dragging

def run_snake_at(x, y, currentPath=None, useMouse=True):
    """Run the snake with the first point forced at (x, y) using current base_x0."""
    if currentPath is None:
        x0 = base_x0.copy()
    else:
        x0 = currentPath.copy()
    if useMouse:
        x0[0] = [x, y]
    # Optimize (your optimize returns (result_path, paths))
    prevFreeze = snk.cfg.fixed_vertices
    if not useMouse:
        snk.cfg.fixed_vertices = []
    result_path, paths = snk.optimize(x0, verbose=False)
    snk.cfg.fixed_vertices = prevFreeze
    xT = paths[-1]
    return xT, paths

def on_key_press(event):
    if event.key == 'enter':
        global _last_run, _last_xy, draw_pts, line, draw_line
        xT, allPaths = run_snake_at(*_last_xy, currentPath=None, useMouse=False)
        # Build segments for the LineCollection (one polyline per iteration)
        if cfg.closed:
            segs = [np.vstack([p, p[0]]) for p in allPaths[50::10]]   # close each loop
        else:
            segs = [p for p in allPaths[50::10]]                      # open polylines

        line.set_data([], [])
        draw_line.set_data([], [])
        # Optional: color by iteration (earlier = lighter, later = darker)
        n = len(segs)
        for iSnake in range(n):
            # Update the collection with all paths
            paths_coll.set_segments([segs[iSnake]])
            colors = "blue"
            paths_coll.set_color(colors)
            paths_coll.set_alpha(1.0)
            paths_coll.set_linewidth(4.5)

            fig.canvas.draw_idle()
            fig.savefig(f"snake-test-{iSnake}.png")
            # plt.pause(1)
            print(f"Saved {iSnake}/{n}")
        paths_coll.set_linewidth(1.0)
        paths_coll.set_alpha(0.6)



def on_mouse_move(event):
    global _last_run, _last_xy, draw_pts
    if event.inaxes != ax or event.xdata is None or event.ydata is None:
        return
    x, y = float(event.xdata), float(event.ydata)

    if is_drawing:
        # add points as the mouse moves to build the initial curve
        if not draw_pts or (abs(x - draw_pts[-1][0]) + abs(y - draw_pts[-1][1])) > 0.75:
            draw_pts.append((x, y))
            xs, ys = zip(*draw_pts)
            draw_line.set_data(xs, ys)
            fig.canvas.draw_idle()
        return

    # normal hover: throttle recomputes
    t = time.monotonic()
    if (t - _last_run) < DEBOUNCE:
        return

    lx, ly = _last_xy
    if lx is not None and (abs(x - lx) < 1.0 and abs(y - ly) < 1.0):
        return

    _last_run = t
    _last_xy = (x, y)

    xT, allPaths = run_snake_at(x, y, currentPath=None)
    # if cfg.closed:
    #     line.set_data(np.r_[xT[:, 0], xT[0, 0]], np.r_[xT[:, 1], xT[0, 1]])
    # else:
    #     line.set_data(xT[:, 0], xT[:, 1])
    # fig.canvas.draw_idle()

    # Build segments for the LineCollection (one polyline per iteration)
    if cfg.closed:
        segs = [np.vstack([p, p[0]]) for p in allPaths[50::50]]   # close each loop
    else:
        segs = [p for p in allPaths[50::50]]                      # open polylines

    # Update the collection with all paths
    paths_coll.set_segments(segs)

    # Optional: color by iteration (earlier = lighter, later = darker)
    n = len(segs)
    colors = [cmap(i/(n-1 if n>1 else 1)) for i in range(n)]
    paths_coll.set_color(colors)

    # Keep the final path emphasized with the existing Line2D
    if cfg.closed:
        line.set_data(np.r_[xT[:, 0], xT[0, 0]], np.r_[xT[:, 1], xT[0, 1]])
    else:
        line.set_data(xT[:, 0], xT[:, 1])

    fig.canvas.draw_idle()

def on_button_press(event):
    # begin drawing a custom initial curve
    global is_drawing, draw_pts
    if event.inaxes != ax or event.xdata is None or event.ydata is None:
        return
    is_drawing = True
    draw_pts = [(float(event.xdata), float(event.ydata))]
    draw_line.set_data([event.xdata], [event.ydata])
    fig.canvas.draw_idle()

def on_button_release(event):
    # finalize the drawn curve, resample to cfg.n_points, set as new base_x0
    global is_drawing, base_x0, draw_pts
    if not is_drawing:
        return
    is_drawing = False

    if len(draw_pts) >= 2:
        pts = np.array(draw_pts, dtype=np.float32)
        # If closed, close the loop if the endpoints are near
        if cfg.closed:
            if np.linalg.norm(pts[0] - pts[-1]) > 1.0:
                pts = np.vstack([pts, pts[0]])
        # resample to uniform spacing with the same open/closed setting
        base_x0 = respaced_polyline(pts, cfg.closed, cfg.n_points)
        cfg.L_target = np.sum(arc_lengths(base_x0, cfg.closed)) * 1.1
        # update preview line to show the new base curve
        if cfg.closed:
            draw_line.set_data(np.r_[base_x0[:,0], base_x0[0,0]], np.r_[base_x0[:,1], base_x0[0,1]])
        else:
            draw_line.set_data(base_x0[:,0], base_x0[:,1])
    else:
        # not enough points: keep previous base_x0
        draw_line.set_data([], [])

    draw_pts = []
    fig.canvas.draw_idle()

# connect events
cid_move   = fig.canvas.mpl_connect('motion_notify_event', on_mouse_move)
cid_press  = fig.canvas.mpl_connect('button_press_event', on_button_press)
cid_release= fig.canvas.mpl_connect('button_release_event', on_button_release)
cid_key    = fig.canvas.mpl_connect('key_press_event', on_key_press)

plt.show()