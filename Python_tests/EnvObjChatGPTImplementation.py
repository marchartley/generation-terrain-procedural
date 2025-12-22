import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# Utility functions
# ============================================================

def finite_gradient(field, dx):
    """Compute central-difference gradient of scalar field.
    Returns (gy, gx) arrays with same shape as field."""
    gy = np.zeros_like(field)
    gx = np.zeros_like(field)

    gy[1:-1, :] = (field[2:, :] - field[:-2, :]) / (2 * dx)
    gx[:, 1:-1] = (field[:, 2:] - field[:, :-2]) / (2 * dx)

    # Neumann boundary: copy neighbours
    gy[0, :] = gy[1, :]
    gy[-1, :] = gy[-2, :]
    gx[:, 0] = gx[:, 1]
    gx[:, -1] = gx[:, -2]
    return gy, gx


def finite_laplacian(field, dx):
    """5-point Laplacian of scalar field."""
    lap = np.zeros_like(field)
    lap[1:-1, 1:-1] = (
        field[2:, 1:-1] + field[:-2, 1:-1] +
        field[1:-1, 2:] + field[1:-1, :-2] -
        4.0 * field[1:-1, 1:-1]
    ) / (dx * dx)
    # Neumann boundary: copy neighbours
    lap[0, :] = lap[1, :]
    lap[-1, :] = lap[-2, :]
    lap[:, 0] = lap[:, 1]
    lap[:, -1] = lap[:, -2]
    return lap


def normalize01(x, eps=1e-8):
    """Normalize array to [0,1] (robust)."""
    mn = np.min(x)
    mx = np.max(x)
    if mx - mn < eps:
        return np.zeros_like(x)
    return (x - mn) / (mx - mn)


def gaussian_kernel_distance(dist2, sigma):
    """Gaussian of squared distance."""
    return np.exp(-dist2 / (2 * sigma * sigma))


# ============================================================
# Materials and Environment
# ============================================================

class MaterialLayer:
    def __init__(self, name, c, rho=0.0, nu=1.0, D=1.0, mu=0.01):
        """
        c  : (ny, nx) scalar field (concentration)
        rho: gravity sensitivity (∇H term)
        nu : flow sensitivity (W term)
        D  : diffusion coefficient
        mu : decay coefficient
        """
        self.name = name
        self.c = c
        self.rho = rho
        self.nu = nu
        self.D = D
        self.mu = mu


class Environment:
    def __init__(self, H, L, W, materials, dx=1.0):
        """
        H : (ny, nx) signed height above water (h - L)
        L : water level (scalar)
        W : (ny, nx, 2) water current vector field
        materials : dict[name -> MaterialLayer]
        dx : grid spacing
        """
        self.H = H
        self.L = L
        self.W = W
        self.materials = materials
        self.dx = dx

    @property
    def shape(self):
        return self.H.shape


# ============================================================
# Environmental Objects
# ============================================================

class EnvObjectClass:
    def __init__(self,
                 name,
                 primitive_type,
                 fitness_fn,
                 skeleton_fit_fn,
                 height_modifier_fn,
                 flow_modifier_fn,
                 material_modifier_fn,
                 fitness_threshold=0.5):
        self.name = name
        self.primitive_type = primitive_type  # 'point' | 'curve' | 'region'
        self.fitness_fn = fitness_fn
        self.skeleton_fit_fn = skeleton_fit_fn
        self.height_modifier_fn = height_modifier_fn
        self.flow_modifier_fn = flow_modifier_fn
        self.material_modifier_fn = material_modifier_fn
        self.fitness_threshold = fitness_threshold


class EnvObjectInstance:
    def __init__(self, obj_class, skeleton, params=None):
        """
        skeleton:
          - point: (y, x)
          - curve: [(y,x), ...]
          - region: mask (ny, nx) boolean
        """
        self.obj_class = obj_class
        self.skeleton = skeleton
        self.params = params or {}


# ============================================================
# Helper: sampling & hill-climbing for points
# ============================================================

def random_candidate_positions(shape, n_samples):
    ny, nx = shape
    ys = np.random.randint(0, ny, size=n_samples)
    xs = np.random.randint(0, nx, size=n_samples)
    return list(zip(ys, xs))


def hill_climb_point(env, fitness_fn, start_pos, max_steps=20, step_radius=1):
    """
    Simple discrete hill-climbing on ω.
    """
    ny, nx = env.shape
    y, x = start_pos
    best_val = fitness_fn(env, y, x)

    for _ in range(max_steps):
        improved = False
        for dy in range(-step_radius, step_radius + 1):
            for dx in range(-step_radius, step_radius + 1):
                yy = y + dy
                xx = x + dx
                if yy < 0 or yy >= ny or xx < 0 or xx >= nx:
                    continue
                val = fitness_fn(env, yy, xx)
                if val > best_val:
                    best_val = val
                    y, x = yy, xx
                    improved = True
        if not improved:
            break
    return (y, x), best_val


# ============================================================
# Skeleton fitting for curves and regions (simplified)
# ============================================================

def trace_downhill(env, start_pos, max_length=200, step_size=1.0):
    """
    Trace a curve along steepest descent of H (for rivers).
    Returns list of (y,x) points.
    """
    H = env.H
    dx = env.dx
    ny, nx = H.shape
    curve = []

    y, x = start_pos
    for _ in range(max_length):
        curve.append((int(round(y)), int(round(x))))
        # Gradient of H
        gy, gx = finite_gradient(H, dx)
        # Negative gradient (downhill)
        vy = -gy[int(round(y)), int(round(x))]
        vx = -gx[int(round(y)), int(round(x))]
        length = np.hypot(vy, vx)
        if length < 1e-5:
            break
        vy /= length
        vx /= length
        y += vy * step_size
        x += vx * step_size
        if y < 0 or y >= ny or x < 0 or x >= nx:
            break
    return curve


def ridge_follow(env, fitness_fn, start_pos, max_length=200, step_size=1.0):
    """
    Simple 'ridge following' snake for reefs & canyons: follow gradient of fitness.
    """
    ny, nx = env.shape
    # Build fitness field locally by sampling; for simplicity, recompute as needed
    y, x = start_pos
    curve = []
    for _ in range(max_length):
        yy = int(round(y))
        xx = int(round(x))
        if yy < 0 or yy >= ny or xx < 0 or xx >= nx:
            break
        curve.append((yy, xx))
        # Local numerical gradient of fitness
        # Use finite differences on 3x3 neighbourhood
        vals = np.zeros((3, 3))
        for j in range(-1, 2):
            for i in range(-1, 2):
                yj = np.clip(yy + j, 0, ny - 1)
                xi = np.clip(xx + i, 0, nx - 1)
                vals[j + 1, i + 1] = fitness_fn(env, yj, xi)
        # gradient in local coordinates (rows= y, cols=x)
        gy_local = (vals[2, 1] - vals[0, 1]) / 2.0
        gx_local = (vals[1, 2] - vals[1, 0]) / 2.0
        length = np.hypot(gy_local, gx_local)
        if length < 1e-5:
            break
        gy_step = gy_local / length
        gx_step = gx_local / length
        y += gy_step * step_size
        x += gx_step * step_size
    return curve


def grow_region_threshold(env, fitness_fn, seed_pos, threshold=0.5, max_size=5000):
    """
    Simple region-growing: BFS from seed while ω >= threshold.
    Returns boolean mask.
    """
    ny, nx = env.shape
    mask = np.zeros((ny, nx), dtype=bool)
    queue = [seed_pos]
    count = 0

    while queue and count < max_size:
        y, x = queue.pop()
        if mask[y, x]:
            continue
        if fitness_fn(env, y, x) < threshold:
            continue
        mask[y, x] = True
        count += 1
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                if dy == 0 and dx == 0:
                    continue
                yy = y + dy
                xx = x + dx
                if 0 <= yy < ny and 0 <= xx < nx and not mask[yy, xx]:
                    queue.append((yy, xx))
    return mask


# ============================================================
# Environmental Modifiers
# ============================================================

def height_modifier_cone(env, center, radius, height_scale):
    """
    Simple conical bump (positive or negative).
    center: (y,x), radius in grid units, height_scale: peak magnitude.
    """
    ny, nx = env.shape
    Y, X = np.indices((ny, nx))
    cy, cx = center
    dist2 = (Y - cy)**2 + (X - cx)**2
    r2 = radius * radius
    mask = dist2 <= r2
    # Cone: h(r) = (1 - r/R)
    r = np.sqrt(dist2[mask])
    h = height_scale * (1.0 - r / radius)
    Hplus = np.zeros_like(env.H)
    Hplus[mask] = h
    return Hplus


def flow_modifier_obstacle(env, skeleton, radius, slowdown_factor=0.5):
    """
    Simple modifier: slow flow magnitude inside radius around skeleton.
    skeleton can be point or list of points.
    """
    ny, nx = env.shape
    Wplus = np.zeros_like(env.W)
    Y, X = np.indices((ny, nx))

    if isinstance(skeleton, tuple):
        points = [skeleton]
    else:
        points = skeleton

    for (cy, cx) in points:
        dist2 = (Y - cy)**2 + (X - cx)**2
        mask = dist2 <= radius * radius
        # scale velocity: (1 - slowdown_factor) * W
        Wplus[mask] += -(1.0 - slowdown_factor) * env.W[mask]
    return Wplus


def material_modifier_gaussian(env, skeleton, mat_name,
                               deposit=0.0, absorb=0.0,
                               sigma=3.0):
    """
    Create source term S for a single material from an object.
    deposit>0 adds material, absorb>0 removes.
    For simplicity: S = deposit * G - absorb * G * c
    Here we only implement deposit; absorption can be done by
    negative deposit or post-processing.
    """
    ny, nx = env.shape
    Y, X = np.indices((ny, nx))
    if isinstance(skeleton, tuple):
        points = [skeleton]
    elif isinstance(skeleton, np.ndarray) and skeleton.dtype == bool:
        # region: approximate as all pixels in region
        ys, xs = np.where(skeleton)
        points = list(zip(ys, xs))
    else:
        # curve: list
        points = skeleton

    S = np.zeros((ny, nx), dtype=np.float32)
    for (cy, cx) in points:
        dist2 = (Y - cy)**2 + (X - cx)**2
        G = gaussian_kernel_distance(dist2, sigma)
        S += deposit * G
    return mat_name, S


# ============================================================
# Advection–Diffusion–Decay solver for materials
# ============================================================

def update_materials(env, objects, dt=0.1, max_iters=20, tol=1e-3):
    """
    For each material, perform a few iterations of explicit
    advection–diffusion–decay with sources from objects.
    This is a simplified steady-state relaxation.
    """
    H = env.H
    W = env.W
    dx = env.dx

    # Precompute gradients of H once
    gHy, gHx = finite_gradient(H, dx)

    # Build source term per material from all objects
    sources = {name: np.zeros_like(layer.c) for name, layer in env.materials.items()}
    for o in objects:
        mat_src = o.obj_class.material_modifier_fn(env, o)
        if mat_src is None:
            continue
        # mat_src can be dict or (mat_name, S)
        if isinstance(mat_src, dict):
            for name, S in mat_src.items():
                sources[name] += S
        else:
            name, S = mat_src
            sources[name] += S

    # Relaxation iterations
    for _ in range(max_iters):
        max_change = 0.0
        for name, layer in env.materials.items():
            c = layer.c
            # velocity for this material
            Phi_y = layer.rho * gHy + layer.nu * W[..., 0]
            Phi_x = layer.rho * gHx + layer.nu * W[..., 1]

            # Upwind-style advection (very crude)
            # Compute divergence of Phi * c ≈ ∇·(Phi * c)
            # First compute fluxes
            # Shifted indices for x direction
            c_right = np.roll(c, -1, axis=1)
            c_left = np.roll(c, 1, axis=1)
            vx = Phi_x
            # simple upwind
            flux_x = np.where(vx > 0, vx * c_left, vx * c_right)

            c_down = np.roll(c, -1, axis=0)
            c_up = np.roll(c, 1, axis=0)
            vy = Phi_y
            flux_y = np.where(vy > 0, vy * c_up, vy * c_down)

            div = (flux_x - np.roll(flux_x, 1, axis=1) +
                   flux_y - np.roll(flux_y, 1, axis=0)) / dx

            lap = finite_laplacian(c, dx)

            S = sources[name]
            # explicit Euler
            dc_dt = -div + layer.D * lap - layer.mu * c + S
            c_new = c + dt * dc_dt
            c_new = np.maximum(c_new, 0.0)  # no negative concentration
            change = np.max(np.abs(c_new - c))
            max_change = max(max_change, change)
            layer.c = c_new
        if max_change < tol:
            break


# ============================================================
# Fitness functions for your 8 object types
# Each fitness takes (env, y, x) and returns scalar in [0,1]
# ============================================================

def depth_at(env, y, x):
    return -env.H[y, x]  # positive below water


def slope_magnitude(env):
    gy, gx = finite_gradient(env.H, env.dx)
    return np.hypot(gy, gx)


def flow_magnitude(env):
    return np.linalg.norm(env.W, axis=-1)


def safe_get_mat(env, name):
    if name in env.materials:
        return env.materials[name].c
    # default zero if missing
    return np.zeros(env.shape)


# -- Coral ---------------------------------------------------

def fitness_coral(env, y, x):
    depth = depth_at(env, y, x)
    flow = flow_magnitude(env)[y, x]
    sand = safe_get_mat(env, "sand")[y, x]
    silt = safe_get_mat(env, "silt")[y, x]
    salt = safe_get_mat(env, "salt")[y, x]
    nutrients = safe_get_mat(env, "nutrients")[y, x]
    # approximate light: high when low silt and shallow
    light = np.exp(-0.1 * depth) * np.exp(-2.0 * silt)

    # Depth preference: 2-25m
    depth_opt = 10.0
    depth_sigma = 8.0
    f_depth = np.exp(-((depth - depth_opt) ** 2) / (2 * depth_sigma ** 2))
    # Flow preference
    flow_opt = 0.5
    flow_sigma = 0.5
    f_flow = np.exp(-((flow - flow_opt) ** 2) / (2 * flow_sigma ** 2))

    # simple combination
    val = f_depth * f_flow * (1 - silt) * salt * light
    # small penalty if too sandy (smothering)
    val *= np.exp(-2.0 * sand)
    # nutrients ok but penalize too high (eutrophication)
    val *= np.exp(-((nutrients - 0.5) ** 2))

    return float(np.clip(val, 0.0, 1.0))


# -- Reef ----------------------------------------------------

def fitness_reef(env, y, x):
    depth = depth_at(env, y, x)
    silt = safe_get_mat(env, "silt")[y, x]
    limestone = safe_get_mat(env, "limestone")[y, x]
    flow = flow_magnitude(env)[y, x]

    depth_opt = 8.0
    depth_sigma = 6.0
    f_depth = np.exp(-((depth - depth_opt) ** 2) / (2 * depth_sigma ** 2))

    flow_opt = 0.7
    flow_sigma = 0.7
    f_flow = np.exp(-((flow - flow_opt) ** 2) / (2 * flow_sigma ** 2))

    val = f_depth * f_flow * (1 - silt) * limestone
    return float(np.clip(val, 0.0, 1.0))


# -- Island --------------------------------------------------

def fitness_island(env, y, x):
    H = env.H[y, x]
    if H <= 0:
        return 0.0  # must be above water
    sand = safe_get_mat(env, "sand")[y, x]
    moisture = safe_get_mat(env, "moisture")[y, x]
    slope = slope_magnitude(env)[y, x]

    f_slope = np.exp(-10.0 * slope)  # penalize steep slopes
    # moisture preference: moderate (~0.5)
    f_moist = np.exp(-((moisture - 0.5) ** 2) / (2 * 0.25 ** 2))
    # print(f"For island: {sand} * {f_slope} * {f_moist}")
    val = sand * f_slope * f_moist
    return float(np.clip(val, 0.0, 1.0))


# -- Algae ---------------------------------------------------

def fitness_algae(env, y, x):
    depth = depth_at(env, y, x)
    nutrients = safe_get_mat(env, "nutrients")[y, x]
    silt = safe_get_mat(env, "silt")[y, x]
    flow = flow_magnitude(env)[y, x]

    light = np.exp(-0.2 * depth) * np.exp(-2.0 * silt)

    depth_opt = 5.0
    depth_sigma = 5.0
    f_depth = np.exp(-((depth - depth_opt) ** 2) / (2 * depth_sigma ** 2))

    flow_opt = 0.3
    flow_sigma = 0.5
    f_flow = np.exp(-((flow - flow_opt) ** 2) / (2 * flow_sigma ** 2))

    val = f_depth * f_flow * light * nutrients
    return float(np.clip(val, 0.0, 1.0))


# -- Lagoon --------------------------------------------------

def fitness_lagoon(env, y, x):
    depth = depth_at(env, y, x)
    flow = flow_magnitude(env)[y, x]
    silt = safe_get_mat(env, "silt")[y, x]
    nutrients = safe_get_mat(env, "nutrients")[y, x]
    salt = safe_get_mat(env, "salt")[y, x]

    # shallow
    depth_opt = 3.0
    depth_sigma = 3.0
    f_depth = np.exp(-((depth - depth_opt) ** 2) / (2 * depth_sigma ** 2))

    # low flow
    f_flow = np.exp(-5.0 * flow)

    # salinity slightly lower than ocean (~0.7-0.9 if ocean ~1.0)
    f_saltmix = np.exp(-((salt - 0.8) ** 2) / (2 * 0.1 ** 2))

    val = f_depth * f_flow * silt * nutrients * f_saltmix
    return float(np.clip(val, 0.0, 1.0))


# -- River ---------------------------------------------------

def fitness_river(env, y, x):
    H = env.H
    slope = slope_magnitude(env)[y, x]
    moisture = safe_get_mat(env, "moisture")[y, x]
    salt = safe_get_mat(env, "salt")[y, x]
    silt = safe_get_mat(env, "silt")[y, x]

    # want fresh / low salinity
    f_salt = np.exp(-10.0 * salt)
    f_moist = np.exp(-((moisture - 0.8) ** 2) / (2 * 0.2 ** 2))
    # prefer some slope
    f_slope = 1.0 - np.exp(-5.0 * slope)

    val = f_salt * f_moist * f_slope * (0.5 + 0.5 * silt)
    return float(np.clip(val, 0.0, 1.0))


# -- Canyon --------------------------------------------------

def fitness_canyon(env, y, x):
    slope = slope_magnitude(env)[y, x]
    flow = flow_magnitude(env)[y, x]
    sand = safe_get_mat(env, "sand")[y, x]
    depth = depth_at(env, y, x)

    # prefer deeper offshore zones
    depth_opt = 20.0
    depth_sigma = 20.0
    f_depth = np.exp(-((depth - depth_opt) ** 2) / (2 * depth_sigma ** 2))

    val = slope * flow * f_depth * np.exp(-2.0 * sand)
    return float(np.clip(val, 0.0, 1.0))


# -- Rock ----------------------------------------------------

def fitness_rock(env, y, x):
    sand = safe_get_mat(env, "sand")[y, x]
    moisture = safe_get_mat(env, "moisture")[y, x]
    slope = slope_magnitude(env)[y, x]

    val = (1 - sand) * (1 - moisture) * (0.5 + 0.5 * slope)
    return float(np.clip(val, 0.0, 1.0))


# ============================================================
# Skeleton fitting wrappers for each class
# ============================================================

# Coral, Rock: point + hill climbing

def skeleton_fit_point(env, fitness_fn, seed_pos):
    pos, _ = hill_climb_point(env, fitness_fn, seed_pos)
    return pos


# Reef & Canyon: curve along gradient of fitness

def skeleton_fit_reef(env, seed_pos):
    return ridge_follow(env, fitness_reef, seed_pos, max_length=150, step_size=1.0)


def skeleton_fit_canyon(env, seed_pos):
    return ridge_follow(env, fitness_canyon, seed_pos, max_length=200, step_size=1.0)


# River: downhill curve

def skeleton_fit_river(env, seed_pos):
    return trace_downhill(env, seed_pos, max_length=250, step_size=1.0)


# Lagoon, Island, Algae: region growing from fitness

def skeleton_fit_region_from_fitness(env, fitness_fn, seed_pos, thr=0.005):
    return grow_region_threshold(env, fitness_fn, seed_pos, threshold=thr, max_size=8000)


# ============================================================
# Per-class modifier functions (very simplified)
# ============================================================

# Height modifiers

def height_mod_coral(env, obj):
    center = obj.skeleton
    return height_modifier_cone(env, center, radius=3, height_scale=0.5)


def height_mod_rock(env, obj):
    center = obj.skeleton
    return height_modifier_cone(env, center, radius=2, height_scale=0.7)


def height_mod_island(env, obj):
    mask = obj.skeleton
    ys, xs = np.where(mask)
    if len(ys) == 0:
        return np.zeros_like(env.H)
    cy = int(np.mean(ys))
    cx = int(np.mean(xs))
    radius = max(5, int(np.sqrt(len(ys)) / 2))
    return height_modifier_cone(env, (cy, cx), radius=radius, height_scale=5.0)


def height_mod_reef(env, obj):
    curve = obj.skeleton
    if not curve:
        return np.zeros_like(env.H)
    # sum of small cones along curve
    Hplus = np.zeros_like(env.H)
    for p in curve[::max(1, len(curve)//20)]:
        Hplus += height_modifier_cone(env, p, radius=2, height_scale=0.5)
    return Hplus


def height_mod_canyon(env, obj):
    curve = obj.skeleton
    if not curve:
        return np.zeros_like(env.H)
    Hplus = np.zeros_like(env.H)
    # negative cones along canyon
    for p in curve[::max(1, len(curve)//20)]:
        Hplus += height_modifier_cone(env, p, radius=3, height_scale=-1.0)
    return Hplus


def height_mod_generic_zero(env, obj):
    return np.zeros_like(env.H)


# Flow modifiers

def flow_mod_obstacle_point(env, obj, radius=4, slowdown=0.5):
    return flow_modifier_obstacle(env, obj.skeleton, radius, slowdown_factor=slowdown)


def flow_mod_obstacle_curve(env, obj, radius=4, slowdown=0.4):
    return flow_modifier_obstacle(env, obj.skeleton, radius, slowdown_factor=slowdown)


def flow_mod_zero(env, obj):
    return np.zeros_like(env.W)


# Material modifiers

def mat_mod_coral(env, obj):
    # Coral deposits limestone
    name, S = material_modifier_gaussian(env, obj.skeleton, "limestone", deposit=0.05, sigma=3.0)
    return {name: S}


def mat_mod_reef(env, obj):
    # Reef deposits limestone, traps silt
    name1, S1 = material_modifier_gaussian(env, obj.skeleton, "limestone", deposit=0.02, sigma=4.0)
    name2, S2 = material_modifier_gaussian(env, obj.skeleton, "silt", deposit=-0.02, sigma=4.0)
    return {name1: S1, name2: S2}


def mat_mod_algae(env, obj):
    # Algae consume nutrients slightly
    name, S = material_modifier_gaussian(env, obj.skeleton, "nutrients", deposit=-0.02, sigma=3.0)
    return {name: S}


def mat_mod_lagoon(env, obj):
    # Lagoon accumulates silt and nutrients
    name1, S1 = material_modifier_gaussian(env, obj.skeleton, "silt", deposit=0.03, sigma=5.0)
    name2, S2 = material_modifier_gaussian(env, obj.skeleton, "nutrients", deposit=0.01, sigma=5.0)
    return {name1: S1, name2: S2}


def mat_mod_river(env, obj):
    # River brings silt and nutrients, reduces salt slightly
    name1, S1 = material_modifier_gaussian(env, obj.skeleton, "silt", deposit=0.03, sigma=3.0)
    name2, S2 = material_modifier_gaussian(env, obj.skeleton, "nutrients", deposit=0.02, sigma=3.0)
    name3, S3 = material_modifier_gaussian(env, obj.skeleton, "salt", deposit=-0.02, sigma=3.0)
    return {name1: S1, name2: S2, name3: S3}


def mat_mod_island(env, obj):
    # Island accumulates sand
    name, S = material_modifier_gaussian(env, obj.skeleton, "sand", deposit=0.05, sigma=6.0)
    return {name: S}


def mat_mod_canyon(env, obj):
    # Canyon reduces sand (erosion)
    name, S = material_modifier_gaussian(env, obj.skeleton, "sand", deposit=-0.03, sigma=4.0)
    return {name: S}


def mat_mod_rock(env, obj):
    # Rocks locally reduce sand
    name, S = material_modifier_gaussian(env, obj.skeleton, "sand", deposit=-0.02, sigma=2.0)
    return {name: S}


def mat_mod_zero(env, obj):
    return None


# ============================================================
# Generator (Algorithm 1-style)
# ============================================================

class EnvObjectGenerator:
    def __init__(self, env, object_classes,
                 max_iters=20,
                 samples_per_class=200):
        self.env = env
        self.object_classes = object_classes
        self.max_iters = max_iters
        self.samples_per_class = samples_per_class
        self.objects = []

    def recompute_height_and_flow(self):
        # base H is whatever env.H currently is; we add modifiers
        H_base = self.env.H.copy()
        H_plus_total = np.zeros_like(H_base)
        W_plus_total = np.zeros_like(self.env.W)
        for o in self.objects:
            H_plus_total += o.obj_class.height_modifier_fn(self.env, o)
            # W_plus_total += o.obj_class.flow_modifier_fn(self.env, o)
        self.env.H = H_base + H_plus_total
        # simple terrain-aware flow: for now, keep original W and add obstacles
        self.env.W = self.env.W + W_plus_total

    def spawn_new_objects(self):
        ny, nx = self.env.shape
        for obj_class in self.object_classes:
            # sample candidate positions
            candidates = random_candidate_positions((ny, nx),
                                                    self.samples_per_class)
            best_pos = None
            best_val = -1.0
            for (y, x) in candidates:
                val = obj_class.fitness_fn(self.env, y, x)
                if val > best_val:
                    best_val = val
                    best_pos = (y, x)
            if best_pos is None or best_val < obj_class.fitness_threshold:
                print(f"Couldn't add {obj_class.name}: best pos = {best_pos}, best val = {best_val} / {obj_class.fitness_threshold}")
                continue

            # skeleton fitting
            if obj_class.primitive_type == "point":
                skel = obj_class.skeleton_fit_fn(self.env, obj_class.fitness_fn, best_pos)
            elif obj_class.primitive_type == "curve":
                # skeleton_fit_fn takes (env, seed_pos)
                skel = obj_class.skeleton_fit_fn(self.env, best_pos)
            elif obj_class.primitive_type == "region":
                # skeleton_fit_fn takes (env, fitness_fn, seed_pos)
                skel = obj_class.skeleton_fit_fn(self.env,
                                                 obj_class.fitness_fn,
                                                 best_pos)
            else:
                raise ValueError("Unknown primitive type")

            new_obj = EnvObjectInstance(obj_class, skel)
            print(f"Added a {obj_class.name}")
            self.objects.append(new_obj)

    def prune_unfit_objects(self):
        kept = []
        for o in self.objects:
            if o.obj_class.primitive_type == "point":
                y, x = o.skeleton
                val = o.obj_class.fitness_fn(self.env, y, x)
            elif o.obj_class.primitive_type == "curve":
                # average fitness along curve
                if not o.skeleton:
                    val = 0.0
                else:
                    vals = [o.obj_class.fitness_fn(self.env, y, x)
                            for (y, x) in o.skeleton]
                    val = float(np.mean(vals))
            else:  # region
                mask = o.skeleton
                ys, xs = np.where(mask)
                if len(ys) == 0:
                    val = 0.0
                else:
                    vals = [o.obj_class.fitness_fn(self.env, int(y), int(x))
                            for y, x in zip(ys, xs)]
                    val = float(np.mean(vals))
            if val >= o.obj_class.fitness_threshold * 0.5:  # some hysteresis
                kept.append(o)
        self.objects = kept

    def run(self):
        for it in range(self.max_iters):
            print(f"Iteration {it+1}/{self.max_iters}, objects: {len(self.objects)}")
            # 1. Instantiate new objects
            self.spawn_new_objects()
            # 2. Recompute environment from objects
            self.recompute_height_and_flow()
            # 3. Update materials with advection–diffusion–decay
            update_materials(self.env, self.objects, dt=0.05, max_iters=10, tol=1e-3)
            # 4. Prune unfit objects
            self.prune_unfit_objects()
        return self.objects


# ============================================================
# Example: minimal setup and run
# ============================================================

def example_demo():
    ny, nx = 128, 128
    dx = 1.0

    # Base height: simple island-ish bump + noise
    Y, X = np.indices((ny, nx))
    cy, cx = ny // 2, nx // 2
    r = np.sqrt((Y - cy)**2 + (X - cx)**2)
    base_h = 5.0 * np.exp(-(r**2) / (2 * (30.0**2)))  # above sea at centre
    base_h += 0.2 * np.random.randn(ny, nx)

    water_level = 2.0
    H = base_h - water_level

    # Simple circular flow around center
    W = np.zeros((ny, nx, 2), dtype=np.float32)
    # tangential
    W[..., 0] = -(X - cx)
    W[..., 1] = (Y - cy)
    norm = np.linalg.norm(W, axis=-1, keepdims=True) + 1e-6
    W = 0.3 * W / norm

    # Materials initialised to simple gradients / noise
    materials = {
        "sand": MaterialLayer("sand", c=normalize01(np.maximum(0, -H)) * 0.5),
        "limestone": MaterialLayer("limestone", c=np.zeros((ny, nx))),
        "silt": MaterialLayer("silt", c=np.zeros((ny, nx)), nu=1.5),
        "nutrients": MaterialLayer("nutrients", c=normalize01(np.random.rand(ny, nx))),
        "salt": MaterialLayer("salt", c=np.ones((ny, nx))),
        "moisture": MaterialLayer("moisture", c=normalize01(np.maximum(0, H))),
    }

    env = Environment(H=H, L=water_level, W=W, materials=materials, dx=dx)

    # Declare object classes
    coral_class = EnvObjectClass(
        name="coral",
        primitive_type="point",
        fitness_fn=fitness_coral,
        skeleton_fit_fn=skeleton_fit_point,
        height_modifier_fn=height_mod_coral,
        flow_modifier_fn=flow_mod_obstacle_point,
        material_modifier_fn=mat_mod_coral,
        fitness_threshold=0.4,
    )

    reef_class = EnvObjectClass(
        name="reef",
        primitive_type="curve",
        fitness_fn=fitness_reef,
        skeleton_fit_fn=skeleton_fit_reef,
        height_modifier_fn=height_mod_reef,
        flow_modifier_fn=flow_mod_obstacle_curve,
        material_modifier_fn=mat_mod_reef,
        fitness_threshold=0.3,
    )

    lagoon_class = EnvObjectClass(
        name="lagoon",
        primitive_type="region",
        fitness_fn=fitness_lagoon,
        skeleton_fit_fn=skeleton_fit_region_from_fitness,
        height_modifier_fn=height_mod_generic_zero,
        flow_modifier_fn=flow_mod_zero,
        material_modifier_fn=mat_mod_lagoon,
        fitness_threshold=0.3,
    )

    island_class = EnvObjectClass(
        name="island",
        primitive_type="region",
        fitness_fn=fitness_island,
        skeleton_fit_fn=skeleton_fit_region_from_fitness,
        height_modifier_fn=height_mod_island,
        flow_modifier_fn=height_mod_generic_zero,
        material_modifier_fn=mat_mod_island,
        fitness_threshold=0.001,
    )

    river_class = EnvObjectClass(
        name="river",
        primitive_type="curve",
        fitness_fn=fitness_river,
        skeleton_fit_fn=skeleton_fit_river,
        height_modifier_fn=height_mod_generic_zero,
        flow_modifier_fn=flow_mod_zero,
        material_modifier_fn=mat_mod_river,
        fitness_threshold=0.3,
    )

    algae_class = EnvObjectClass(
        name="algae",
        primitive_type="region",
        fitness_fn=fitness_algae,
        skeleton_fit_fn=skeleton_fit_region_from_fitness,
        height_modifier_fn=height_mod_generic_zero,
        flow_modifier_fn=flow_mod_zero,
        material_modifier_fn=mat_mod_algae,
        fitness_threshold=0.3,
    )

    canyon_class = EnvObjectClass(
        name="canyon",
        primitive_type="curve",
        fitness_fn=fitness_canyon,
        skeleton_fit_fn=skeleton_fit_canyon,
        height_modifier_fn=height_mod_canyon,
        flow_modifier_fn=flow_mod_obstacle_curve,
        material_modifier_fn=mat_mod_canyon,
        fitness_threshold=0.3,
    )

    rock_class = EnvObjectClass(
        name="rock",
        primitive_type="point",
        fitness_fn=fitness_rock,
        skeleton_fit_fn=skeleton_fit_point,
        height_modifier_fn=height_mod_rock,
        flow_modifier_fn=flow_mod_obstacle_point,
        material_modifier_fn=mat_mod_rock,
        fitness_threshold=0.4,
    )

    obj_classes = [
        coral_class,
        reef_class,
        lagoon_class,
        island_class,
        river_class,
        algae_class,
        canyon_class,
        rock_class,
    ]

    gen = EnvObjectGenerator(env, obj_classes,
                             max_iters=10,
                             samples_per_class=100)
    objects = gen.run()
    print("Generation finished, number of objects:", len(objects))

    # Optionally: return env and objects for visualization
    return env, objects

# ============================================================
# Matplotlib viewer
# ============================================================

def plot_field(ax, data, title):
    """Display a scalar field with imshow."""
    im = ax.imshow(data, origin="lower", interpolation="nearest")
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
def overlay_objects(ax, objects):
    """
    Overlay environmental objects on a given axes.
    - Points: markers
    - Curves: lines
    - Regions: contours
    """
    for o in objects:
        name = o.obj_class.name.lower()
        if o.obj_class.primitive_type == "point":
            y, x = o.skeleton
            ax.plot(x, y, "o", markersize=4, label=name)
        elif o.obj_class.primitive_type == "curve":
            if not o.skeleton:
                continue
            ys = [p[0] for p in o.skeleton]
            xs = [p[1] for p in o.skeleton]
            ax.plot(xs, ys, "-", linewidth=1.0, label=name)
        elif o.obj_class.primitive_type == "region":
            mask = o.skeleton
            if mask is None or not mask.any():
                continue
            # contour of the region mask
            ax.contour(mask.astype(float), levels=[0.5], origin="lower", linewidths=1.0)

    # Avoid cluttering legend: keep unique labels
    handles, labels = ax.get_legend_handles_labels()
    unique = {}
    for h, l in zip(handles, labels):
        if l not in unique:
            unique[l] = h
    if unique:
        ax.legend(unique.values(), unique.keys(), fontsize="small", loc="upper right")
def view_scene(env, objects):
    """
    Simple overview viewer:
    - Height field
    - Sand
    - Limestone
    - Nutrients
    with objects overlaid on the height field.
    """
    H = env.H
    sand = env.materials.get("sand").c if "sand" in env.materials else None
    limestone = env.materials.get("limestone").c if "limestone" in env.materials else None
    nutrients = env.materials.get("nutrients").c if "nutrients" in env.materials else None
    silt = env.materials.get("silt").c if "silt" in env.materials else None

    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    axH, axSand, axLime, axNutr, axSilt, axEmpty = axes.ravel()

    # Height
    plot_field(axH, H, "Height (H)")
    overlay_objects(axH, objects)

    # Sand
    if sand is not None:
        plot_field(axSand, sand, "Sand")
    else:
        axSand.set_title("Sand (none)")
        axSand.axis("off")

    # Limestone
    if limestone is not None:
        plot_field(axLime, limestone, "Limestone")
    else:
        axLime.set_title("Limestone (none)")
        axLime.axis("off")

    # Nutrients
    if nutrients is not None:
        plot_field(axNutr, nutrients, "Nutrients")
    else:
        axNutr.set_title("Nutrients (none)")
        axNutr.axis("off")

    # Silt
    if silt is not None:
        plot_field(axSilt, silt, "Silt")
    else:
        axSilt.set_title("Silt (none)")
        axSilt.axis("off")

    # Last panel: empty or something else (e.g. moisture or flow magnitude)
    moisture = env.materials.get("moisture").c if "moisture" in env.materials else None
    if moisture is not None:
        plot_field(axEmpty, moisture, "Moisture")
    else:
        axEmpty.set_title("Moisture (none)")
        axEmpty.axis("off")

    fig.tight_layout()
    return fig

if __name__ == "__main__":
    env, objects = example_demo()
    fig = view_scene(env, objects)
    plt.show()
