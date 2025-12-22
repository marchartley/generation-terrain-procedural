import functools
import os
import math
import random
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from progressbar import GranularBar, AnimatedMarker, CurrentTime
from scipy.ndimage import gaussian_filter
import progressbar
import matplotlib.image

originalAx: matplotlib.image.AxesImage = None
erodedAx: matplotlib.image.AxesImage = None
differenceAx: matplotlib.image.AxesImage = None


_worker_obj = None

def _init_worker(obj):
    # Runs once in each process: stash your object (or the parts it needs)
    global _worker_obj
    _worker_obj = obj

def _run_particle(args):
    # args: (x, y, z)
    x, y, z = args
    p = np.array([x, y, z], dtype=np.float32)
    v = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    C = 0.0
    # Expecting:  _, _, _, poses, changes
    return _worker_obj.step_particle(p, v, C)

# -----------------------------
# Utilities
# -----------------------------
def load_heightmap(path_candidates):
    for p in path_candidates:
        if os.path.exists(p):
            img = Image.open(p).convert('L')
            arr = np.asarray(img, dtype=np.float32) / 255.0  # normalize 0..1
            arr *= 100.0
            return arr, p
    # Fallback: synthesize a hill-valley terrain so the demo always runs
    H, W = 256, 256
    y, x = np.mgrid[0:H, 0:W]
    x = (x - W/2) / (W/2)
    y = (y - H/2) / (H/2)
    r = np.sqrt(x*x + y*y)
    z = 0.5*np.exp(-3*r*r) + 0.25*np.exp(-20*((x+0.2)**2 + (y-0.1)**2))
    z = (z - z.min()) / (z.max() - z.min())
    return z.astype(np.float32), None

def save_heightmap(arr, path_png, path_npy=None, diff=None):
    arr_clipped = np.clip(arr, 0.0, 1.0)
    img = Image.fromarray((arr_clipped * 255.0).astype(np.uint8), mode='L')
    img.save(path_png)
    # np.save(path_npy, arr_clipped)
    if diff is not None:
        # Create a TwoSlopeNorm centered at zero
        norm = matplotlib.colors.TwoSlopeNorm(vmin=np.min(diff) - 0.01, vcenter=0, vmax=np.max(diff) + 0.01)
        cmap = plt.get_cmap('coolwarm')
        diff_rgb = cmap(norm(diff))[:, :, :3]
        img_rgb = Image.fromarray((diff_rgb * 255).astype(np.uint8), mode='RGB')
        img_rgb.save(path_png.replace('.png', '_diff.png'))

# Bilinear sample of height field
def sample_h(h, x, y):
    H, W = h.shape
    if x < 0 or y < 0 or x > W-1 or y > H-1:
        return None
    x0 = int(np.floor(x)); x1 = min(x0+1, W-1)
    y0 = int(np.floor(y)); y1 = min(y0+1, H-1)
    tx = x - x0; ty = y - y0
    z00 = h[y0, x0]; z10 = h[y0, x1]
    z01 = h[y1, x0]; z11 = h[y1, x1]
    z0 = z00*(1-tx) + z10*tx
    z1 = z01*(1-tx) + z11*tx
    return z0*(1-ty) + z1*ty

# Central differences normal (unnormalized gradient)
def surface_normal(h, x, y):
    H, W = h.shape
    x0 = int(np.clip(round(x), 1, W-2))
    y0 = int(np.clip(round(y), 1, H-2))
    # dzdx = (h[y0, x0+1] - h[y0, x0-1]) * 0.5
    # dzdy = (h[y0+1, x0] - h[y0-1, x0]) * 0.5
    dzdx = (h[y0, x0+1] - h[y0, x0])
    dzdy = (h[y0+1, x0] - h[y0, x0])
    # Height field surface normal for z=h(x,y): n ∝ (-dz/dx, -dz/dy, 1)
    n = np.array([-dzdx, -dzdy, 1.0], dtype=np.float32)
    n /= (np.linalg.norm(n) + 1e-8)
    return n

# Apply half-sphere brush on the height field at contact (cx,cy), radius r_px, with net volume DeltaV (in "height*px^2" units)
def apply_half_sphere_brush(h: np.ndarray, c: np.ndarray, r_px: float, v: np.ndarray, DeltaV: float):
    H, W = h.shape
    cx, cy = c[0], c[1]
    r = int(max(1, round(r_px)))
    x0 = int(max(0, math.floor(cx - r)))
    x1 = int(min(W-1, math.ceil(cx + r)))
    y0 = int(max(0, math.floor(cy - r)))
    y1 = int(min(H-1, math.ceil(cy + r)))
    if x0 >= x1 or y0 >= y1:
        return
    # Precompute normalization: hemisphere volume in pixel-units
    V_hemi = (2.0/3.0) * math.pi * (r_px**3)
    if V_hemi <= 0:
        return
    scale = (DeltaV / V_hemi)
    xx = np.arange(x0, x1+1, dtype=np.float32)
    yy = np.arange(y0, y1+1, dtype=np.float32)
    X, Y = np.meshgrid(xx, yy)
    D2 = (X - cx)**2 + (Y - cy)**2
    mask = D2 <= (r_px**2)
    # sqrt term for half-sphere cap
    delta = np.zeros_like(D2, dtype=np.float32)
    delta[mask] = scale * np.sqrt(np.maximum(r_px**2 - D2[mask], 0.0))
    # Update height field
    h[y0:y1+1, x0:x1+1] += delta

def apply_gaussian_brush(h: np.ndarray, c: np.ndarray, r_px: float, v: np.ndarray, DeltaV: float):
    H, W = h.shape
    cx, cy = c[0], c[1]
    r = int(max(1, round(r_px)))
    x0 = int(max(0, math.floor(cx - r)))
    x1 = int(min(W-1, math.ceil(cx + r)))
    y0 = int(max(0, math.floor(cy - r)))
    y1 = int(min(H-1, math.ceil(cy + r)))
    if x0 >= x1 or y0 >= y1:
        return
    sigma = 0.5 * float(r_px)
    if sigma <= 0.0:
        return

    xx = np.arange(x0, x1 + 1, dtype=np.float32)
    yy = np.arange(y0, y1 + 1, dtype=np.float32)
    X, Y = np.meshgrid(xx, yy)
    D2 = (X - cx) ** 2 + (Y - cy) ** 2

    # Truncate to radius r_px
    mask = D2 <= (r_px ** 2)
    weights = np.zeros_like(D2, dtype=np.float32)
    weights[mask] = np.exp(-0.5 * D2[mask] / (sigma * sigma))

    S = float(weights.sum())
    if S <= 1e-12:
        return

    # Scale so that the discrete integral (sum over pixels) matches DeltaV
    scale = DeltaV / S
    h[y0:y1 + 1, x0:x1 + 1] += scale * weights

def apply_half_ellipsoid_brush(h: np.ndarray, c: np.ndarray, r_px: float,
                               v: np.ndarray, DeltaV: float, k: float = 2.0):
    H, W = h.shape
    cx, cy = float(c[0]), float(c[1])

    # Semi-axes: elongated along v by factor k
    a = r_px * k         # along v
    b = r_px / k         # perpendicular to v
    c_vert = r_px        # vertical semi-axis (keeps baseline volume same as hemisphere)

    # Bounding box (AABB of rotated ellipse)
    rad = int(max(1, round(max(a, b))))
    x0 = max(0, math.floor(cx - rad)); x1 = min(W-1, math.ceil(cx + rad))
    y0 = max(0, math.floor(cy - rad)); y1 = min(H-1, math.ceil(cy + rad))
    if x0 >= x1 or y0 >= y1: return

    # Orientation from v
    vx, vy = float(v[0]), float(v[1])
    n = math.hypot(vx, vy)
    if n < 1e-8:
        ux, uy = 1.0, 0.0   # default orientation
    else:
        ux, uy = vx / n, vy / n
    wx, wy = -uy, ux

    # Volume normalization (same as hemisphere if a*b*c_vert = r^3)
    V_half = (2.0/3.0) * math.pi * (a * b * c_vert)
    if V_half <= 0: return
    scale = DeltaV / V_half

    # Grid
    xx = np.arange(x0, x1+1, dtype=np.float32)
    yy = np.arange(y0, y1+1, dtype=np.float32)
    X, Y = np.meshgrid(xx, yy)

    # Local rotated coordinates
    DX = X - cx; DY = Y - cy
    Xp = DX * ux + DY * uy     # along v (major axis)
    Yp = DX * wx + DY * wy     # perpendicular (minor axis)

    # Elliptic mask and height
    q = (Xp / a)**2 + (Yp / b)**2
    mask = q <= 1.0

    delta = np.zeros_like(q, dtype=np.float32)
    # z = c_vert * sqrt(1 - q), then scale for target DeltaV
    delta[mask] = scale * c_vert * np.sqrt(np.maximum(1.0 - q[mask], 0.0))

    # print(f"Diff = {np.sum(delta)}")

    h[y0:y1+1, x0:x1+1] += delta


# -----------------------------
# Erosion simulation (height field + constant flow)
# -----------------------------
class ErosionSim2p5D:
    def __init__(self, heightmap, rng=np.random.default_rng(42)):
        self.h = heightmap
        self.H, self.W = self.h.shape
        self.rng = rng

        # --- Physical-ish parameters (dimensionless, tuned for demo) ---
        self.g = np.array([0.0, 0.0, -2.0], dtype=np.float32)  # gravity down in z (dimensionless scale)
        self.mu = 0.08      # viscosity proxy
        self.rho_f = 1.0    # medium density (air-like)
        self.flow_u = np.array([0.0, -10.0, 0.0], dtype=np.float32)  # constant flow in +x

        # --- Particle parameters (defaults; can be overridden per scenario) ---
        self.radius_px = 20.0
        self.rho_p = 1000.0     # particle density (heavy)
        self.e_restitution = 0.9
        self.C_max = 100.0 #0.3 * (2.0/3.0)*math.pi*(self.radius_px**3) * 0.5  # capacity ~ half hemisphere volume
        self.tau_c = 1.00      # critical shear stress
        self.n_shear = 0.50     # shear-thinning index
        self.kappa_erosion = 0.6  # erosion strength
        self.lambda_deposit = 0.8  # deposition factor
        self.erosion_scale = 2.0
        self.alpha_flow = 1.0  # scale for flow
        self.amax = 10.0
        self.vmax = 10.0

        self.water_level = 0.0


        # Numerics
        self.dt = 0.1
        self.max_steps = 1000    # steps per particle
        self.contact_tol = 1e-8

        self.simu_acceleration = 1.0

        self.scale_z = 1.#/100.

    def step_particle(self, p, v, C):
        """Integrate one particle for max_steps or until termination.
           Returns updated heightmap via in-place edits and final state.
        """
        last_collision = np.array([np.inf, np.inf, np.inf])
        h0 = self.h.copy()
        newH = self.h.copy()
        positions = [p.copy()]
        center = np.array([h0.shape[0], h0.shape[1], 0]) / 2
        radius = self.radius_px + random.randint(-15, 5)
        # fluid_density = self.rho_f
        for _ in range(self.max_steps):

            u = self.alpha_flow * self.flow_u

            if p[2] > self.water_level:
                fluid_density = 1.0
            else:
                fluid_density = 1000.0
                toCenter = (center - p * np.array([1, 1, 0]))
                u += 2.0 * toCenter / np.linalg.norm(toCenter) #self.flow_u
            # Flow sample (constant here)

            # Body force with buoyancy; use hindrance factor on slip (capacity effect)
            hind = 1.0 # max(0.0, 1.0 - C / max(self.C_max, 1e-8)) ** 5.0  # Richardson–Zaki n≈5
            F_body = hind * (self.rho_p - fluid_density) * self.g

            # Stokes drag towards flow (no vertical target)
            F_drag = 6.0*math.pi*self.mu*radius * (u - v)

            # Velocity-Verlet style semi-implicit update
            a = F_body + F_drag
            if np.linalg.norm(a) > self.amax:
                a = (a * self.amax) / (np.linalg.norm(a))
            v = v + self.dt * a
            v[2] *= self.scale_z
            if np.linalg.norm(v) > self.vmax:
                v = (v * self.vmax) / (np.linalg.norm(v))
            p = p + self.dt * v * self.simu_acceleration
            # print(np.linalg.norm(v))

            positions += [p.copy()]

            # Bounds check (xy within domain, z reasonable)
            if p[0] < 0 or p[0] > self.W-1 or p[1] < 0 or p[1] > self.H-1:
                break
            zsurf = sample_h(h0, p[0], p[1])
            if zsurf is None:
                break
            # Collision with height field surface z = h(x,y)
            if p[2] <= zsurf + self.contact_tol:
                # Project to surface
                p[2] = zsurf
                n = surface_normal(h0, p[0], p[1])

                # Relative tangential speed magnitude for shear proxy
                v_rel = v - np.array([u[0], u[1], 0.0], dtype=np.float32)
                # Remove normal component for shear (tangent slip)
                v_rel_t = v_rel - np.dot(v_rel, n)*n
                gamma_dot = np.linalg.norm(v_rel_t)  # l ~ 1 px
                tau = (gamma_dot**self.n_shear)  # K≈1 absorbed into kappa

                # Detachment rate (only if above threshold)
                E_rate = self.kappa_erosion * max(tau - self.tau_c, 0.0)
                # Hemisphere volume scale for 2.5D
                V_hemi = (2.0/3.0)*math.pi*(radius**3)
                E = E_rate * V_hemi
                # Clamp by remaining capacity
                E = max(0.0, min(E, self.C_max - C))

                # Deposition heuristic using Stokes terminal speed (as weight only)
                ws = (2.0/9.0) * abs(self.g[2]) * (radius**2) * ((self.rho_p - self.rho_f) / max(self.mu,1e-6)) * hind
                D = self.lambda_deposit * abs(ws)
                D = max(0.0, min(D, C))

                # Collision response: restitution on normal component, damp tangential a bit
                vn = np.dot(v, n) * n
                vt = v - vn

                if np.isinf(last_collision[0]):
                    last_collision = p.copy()
                    # last_collision = p.copy() - np.array([1000, 1000, 1000])

                if np.linalg.norm(last_collision - p) > radius * 0.1:
                    DeltaV = D - E  # net terrain change (height*px^2 units)
                    # DeltaV =  min(abs(DeltaV), 10.0) * (-1.0 if DeltaV < 0 else 1.0)
                    # apply_half_sphere_brush(newH, p, radius, v, DeltaV * self.erosion_scale)
                    apply_half_ellipsoid_brush(newH, p, radius, v, DeltaV * self.erosion_scale, 1.0)
                    # apply_half_ellipsoid_brush(newH, p, radius, v, DeltaV * self.erosion_scale, 1.0 + np.clip(np.linalg.norm(v) / 5.0, 0, radius))
                    # print(np.linalg.norm(vt))

                    # Update carried load
                    C = C + E - D
                    last_collision = p.copy()

                # v = (-self.e_restitution) * vn #+ 0.8 * vt  # simple tangential damping
                v = (v - 2 * np.dot(v, n) * n) * self.e_restitution
                # Small lift to avoid re-sticking
                p[2] += 1e-2

            # Termination if almost stopped & high above/below
            # if np.linalg.norm(v) < 1e-5:
            #     break

        if self.lambda_deposit > 0 and not (p[0] < 0 or p[0] > self.W-1 or p[1] < 0 or p[1] > self.H-1) and abs(p[2] - sample_h(h0, p[0], p[1])) < 0.10:
            apply_gaussian_brush(newH, p, radius, np.array([0, 0, 0]), C * self.erosion_scale)
        # print(C)
        return p, v, C, positions, newH - self.h

    def run(self, particles_per_iter=200, iterations=10, initial_map:np.ndarray = None, spawn_mode="rain"):
        H, W = self.H, self.W
        rng = self.rng
        t0 = time.time()
        positions = [[] for _ in range(iterations)]

        bar = progressbar.ProgressBar(max_value=iterations)

        for it in range(iterations):
            bar.update(it)

            # if it % 100 == 0:
            #     save_heightmap(self.h / 100.0, f"/media/marc/Data/NN Datasets/1/result_height--it-{it}.png", diff = self.h - initial_map)

            # Simple "rain": drop particles above random (x,y)
            xs = rng.uniform(0, (W-1)/2, size=particles_per_iter).astype(np.float32)
            ys = rng.uniform(0, H-1, size=particles_per_iter).astype(np.float32)
            zs = np.array([sample_h(self.h, x, y) for x,y in zip(xs,ys)], dtype=np.float32) + (self.contact_tol * 2) #rng.uniform(0.5, 1.0, size=particles_per_iter).astype(np.float32)
            positions[it] = [[] for _ in range(particles_per_iter)]
            # allChanges = np.zeros((*self.h.shape, particles_per_iter))
            combinedChanges = np.zeros_like(self.h)

            with ProcessPoolExecutor(
                    initializer=_init_worker, initargs=(self,), max_workers=None
            ) as ex:
                futures = [ex.submit(_run_particle, (x, y, z))
                           for x, y, z in zip(xs, ys, zs)]

                for i, fut in enumerate(as_completed(futures)):
                    _, _, _, poses, changes = fut.result()
                    positions[it][i] = poses
                    combinedChanges += changes.astype(np.float32)

            # for i in range(particles_per_iter):
            #     _init_worker(self)
            #     p = np.array([xs[i], ys[i], zs[i]], dtype=np.float32)
            #     _, _, _, poses, changes = _run_particle([p[0], p[1], p[2]])
            #     # v = np.array([0.0, 0.0, 0.0], dtype=np.float32)
            #     # C = 0.0  # carried sediment starts empty
            #     # _, _, _, poses, changes = self.step_particle(p, v, C)
            #     positions[it][i] = poses
            #     combinedChanges += changes

            # combinedChanges = np.clip(combinedChanges, -1.0, 1.0)
            combinedChanges = gaussian_filter(combinedChanges, sigma=self.radius_px * 0.25)
            self.h += combinedChanges
            if erodedAx is not None:
                erodedAx.set_data(self.h)
            if differenceAx is not None:
                diff = self.h - initial_map
                norm = matplotlib.colors.TwoSlopeNorm(vmin=np.min(diff) - 0.01, vcenter=0, vmax=np.max(diff) + 0.01)
                differenceAx.set_data(diff)
                differenceAx.set_cmap('bwr')
                differenceAx.set_norm(norm)
                # differenceAx.set_clim(np.min(diff), np.max(diff))
            # plt.pause(0.001)
            plt.suptitle(f"{(it+1)} / {iterations}")
            plt.gcf().canvas.draw_idle()
            plt.gcf().canvas.start_event_loop(0.001)
        bar.finish()
        t1 = time.time()
        return t1 - t0, positions

# -----------------------------
# Main
# -----------------------------
heightmap_paths = [
    "heightmap2.png"
]
h0, used_path = load_heightmap(heightmap_paths)
h_before = h0.copy()

sim = ErosionSim2p5D(h0)


# -----------------------------
# Preview before/after
# -----------------------------
plt.ion()
fig = plt.figure(figsize=(10,4))
plt.subplot(1,3,1)
plt.title("Original heightmap")
originalAx = plt.imshow(h_before, cmap='gray', vmin=0.0, vmax=100.0)
plt.axis('off')

plt.subplot(1,3,2)
plt.title("Eroded heightmap")
erodedAx = plt.imshow(sim.h, cmap='gray', vmin=0.0, vmax=100.0)
plt.axis('off')
plt.tight_layout()

plt.subplot(1,3,3)
plt.title("Difference")
differenceAx = plt.imshow(sim.h - h_before, cmap='gray')
plt.colorbar(label='Value')
plt.axis('off')
plt.tight_layout()
plt.show()
plt.pause(0.01)

elapsed, positions = sim.run(particles_per_iter=20, iterations=1000, initial_map = h_before)

out_png = "/media/marc/Data/NN Datasets/1/result_height.png" #"eroded_heightmap.png"
out_npy = "eroded_heightmap.npy"
# save_heightmap(sim.h / 100.0, out_png, out_npy, diff=sim.h - h_before)


print("Erosion complete in {:.2f}s".format(elapsed))
print("Input used:", used_path if used_path is not None else "(synthetic demo terrain)")
print("Saved outputs:")
print(" -", out_png)
print(" -", out_npy)

plt.ioff()
plt.show()
#
# # --------- Flatten the nested positions list ---------
# def flatten_positions(pos_nested):
#     xs, ys, zs = [], [], []
#     per_particle = []
#     for it in range(len(pos_nested)):
#         row = pos_nested[it]
#         for pid in range(len(row)):
#             # p = np.asarray(row[pid]).reshape(-1)
#             per_particle.append(np.array(row[pid]))
#             for p in row[pid]:
#                 if p.size >= 3 and np.all(np.isfinite(p[:3])):
#                     xs.append(float(p[0])); ys.append(float(p[1])); zs.append(float(p[2]))
#     return np.array(xs), np.array(ys), np.array(zs), per_particle
#
# xs, ys, zs, per_particle = flatten_positions(positions)
#
# # --------- Build a decimated surface for the height field ---------
# H, W = sim.h.shape
# stride = max(1, int(max(H, W)//128))  # adaptively decimate for speed
# ys_grid = np.arange(0, H, stride, dtype=np.int32)
# xs_grid = np.arange(0, W, stride, dtype=np.int32)
# X, Y = np.meshgrid(xs_grid.astype(np.float32), ys_grid.astype(np.float32))
# Z = sim.h[np.ix_(ys_grid, xs_grid)]
#
# # --------- Single 3D plot (no subplots, no explicit colors) ---------
# fig = plt.figure(figsize=(8, 6))
# ax = fig.add_subplot(111, projection='3d')
#
# ax.plot_surface(X, Y, Z, linewidth=0, antialiased=True, alpha=0.7)
# if xs.size > 0:
#     # Scatter is robust even if particles are different per iteration
#     # ax.scatter(xs, ys, zs, s=6, c='red')
#     for particle in per_particle:
#         ax.plot(xs=particle[:, 0], ys=particle[:, 1], zs=particle[:,2])
#
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')
#
# ax.set_xlim(0, W-1)
# ax.set_ylim(0, H-1)
# zmin = float(np.nanmin(Z)); zmax = float(np.nanmax(Z))
# if zs.size > 0:
#     zmin = min(zmin, float(np.nanmin(zs)))
#     zmax = max(zmax, float(np.nanmax(zs)))
# ax.set_zlim(zmin, zmax)
#
# plt.show()