import numpy as np
from scipy.ndimage import maximum_filter, laplace, distance_transform_edt, gaussian_filter
from scipy.optimize import newton
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# --- diffusion + decay ---
def diffuse_decay_step(C, D, lam, dx, dy, dt, bc="neumann"):
    if bc == "neumann":
        Cpad = np.pad(C, 1, mode="edge")
    elif bc == "dirichlet":
        Cpad = np.pad(C, 1, mode="constant", constant_values=0.0)
    else:
        raise ValueError("bc must be 'neumann' or 'dirichlet'")

    lap = ((Cpad[2:, 1:-1] - 2*Cpad[1:-1, 1:-1] + Cpad[:-2, 1:-1]) / dx**2
         + (Cpad[1:-1, 2:] - 2*Cpad[1:-1, 1:-1] + Cpad[1:-1, :-2]) / dy**2)


    return C + dt * (D * lap - lam * C)


def advect_upwind(C, ux, uy, dx, dy, dt, bc="neumann"):
    """
    2D first-order upwind advection with spatially varying ux, uy.
    """
    ny, nx = C.shape

    # --- boundary padding ---
    if bc == "periodic":
        Cx_minus = np.roll(C,  1, axis=1)
        Cx_plus  = np.roll(C, -1, axis=1)
        Cy_minus = np.roll(C,  1, axis=0)
        Cy_plus  = np.roll(C, -1, axis=0)

    elif bc == "neumann":
        Cpad = np.pad(C, ((1, 1), (1, 1)), mode="edge")
        Cx_minus = Cpad[1:-1, 0:-2]
        Cx_plus  = Cpad[1:-1, 2:  ]
        Cy_minus = Cpad[0:-2, 1:-1]
        Cy_plus  = Cpad[2:  , 1:-1]
    else:
        raise ValueError("bc must be 'neumann' or 'periodic'")

    # --- upwind derivative in x ---
    dCdx_pos = (C - Cx_minus) / dx
    dCdx_neg = (Cx_plus - C) / dx
    dCdx = np.where(ux > 0, dCdx_pos, dCdx_neg)

    # --- upwind derivative in y ---
    dCdy_pos = (C - Cy_minus) / dy
    dCdy_neg = (Cy_plus - C) / dy
    dCdy = np.where(uy > 0, dCdy_pos, dCdy_neg)

    # assemble update
    return C - dt * (ux * dCdx + uy * dCdy)


# circular mask
def mask(M, x, y, r):
    ny, nx = M.shape
    Y, X = np.ogrid[:ny, :nx]
    return (X - x)**2 + (Y - y)**2 <= r**2


# --- simulation setup with two sources ---
nx, ny = 128, 128
dx = dy = 1.0
D = 0.9
lam = 0.0
# lam = 0.01
dt = 0.01
t = 100.0
q = 10.0

x = np.arange(nx) * dx
ux = np.ones((ny, nx)) * 1.0        # constant flow to the right
uy = np.zeros((ny, nx))             # default vertical flow = 0
# uy[:, 100:] = 1.0                   # for x > 100, vertical flow downward
# uy_1d = ((x - nx/2)/(nx/3))**3
# make it 2D: same profile in y for each row
# uy = np.broadcast_to(uy_1d, (ny, nx))
# ux_1d = np.sign((nx/2) - x) * -1.0
# ux = np.broadcast_to(ux_1d, (ny, nx))

C = np.zeros((ny, nx), dtype=float)
sources = [mask(C, nx//2, ny//2, 4)]#, mask(C, nx//2, ny//2 + 20, 8)]
sinks   = []#[mask(C, nx//2 + 15, ny//2 + 20, 6)]

nsteps = int(t / dt)

# --- animation setup ---
fig, ax = plt.subplots()
im = ax.imshow(C, origin="lower", vmin=0, vmax=q, animated=True)
cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Concentration")
title = ax.set_title("") # ("ARD with upwind advection, step 0")

# If you want to skip frames to speed up the animation, increase this
steps_per_frame = 100

def step_once(C):
    """Advance the simulation by a single time step."""
    # continuous injection
    for M in sources:
        C[M] = q
    for M in sinks:
        C[M] = 0.0

    # diffusion + decay
    C = diffuse_decay_step(C, D, lam, dx, dy, dt, bc="neumann")

    # advection
    C = advect_upwind(C, ux, uy, dx, dy, dt, bc="neumann")
    return C

def update(frame):
    """Update function for FuncAnimation."""
    global C
    # take several small time steps between frames for smoother physics
    for _ in range(steps_per_frame):
        C = step_once(C)

    im.set_array(C)
    # title.set_text(f"ARD with upwind advection, step {frame * steps_per_frame}")
    return im, title

anim = FuncAnimation(
    fig,
    update,
    frames=nsteps // steps_per_frame,
    interval=30,   # milliseconds between frames (controls playback speed)
    blit=True
)

plt.show()
C[:,:] = 0.0
# If you want to save as mp4 or gif, uncomment one:
# anim.save("ard_process.mp4", fps=30)
anim.save("ard_process.gif", fps=20)

exit()



















import math
import time

import numpy as np
import matplotlib.pyplot as plt

from mpmath import lambertw        # or use scipy.special.lambertw
from scipy.special import k0, k1
from scipy.optimize import newton   # or root_scalar for Brent
#
# gamma = 0.5772156649015328606
# def function(x):
#     return k0(x)
# def functionP(x):
#     return k1(x)
#
# def k0_inv(y, func=k0, funcP=k1, tol=1e-1, maxiter=50000):
#     if y <= 0:
#         raise ValueError("Inverse only defined for y > 0.")
#     # initial guess
#     if y >= 1e0:
#         r0 = 2.0*np.exp(-(y + gamma))
#     elif y <= 1e-6:
#         r0 = 0.5*float(lambertw(np.pi / (y*y)).real)
#     else:
#         # pick one of the asymptotics; both work fine
#         r0 = 0.5*float(lambertw(np.pi / (y*y)).real)
#
#     # Newton step using k1 for derivative
#     f  = lambda r: func(r) - y
#     fp = lambda r: -funcP(r)
#     r  = newton(f, r0, fprime=fp, tol=tol, maxiter=maxiter)
#
#     # Safety: ensure positive root
#     if r <= 0:
#         # fallback to a bracketing method if something odd happens
#         from scipy.optimize import root_scalar
#         # bracket by expanding until sign change
#         a, b = 1e-16, max(r0, 1.0)
#         while func(b) > y:
#             b *= 2.0
#         sol = root_scalar(lambda t: func(t) - y, bracket=(a, b), method='brentq', xtol=tol)
#         r = sol.root
#     return r
#
# x = np.arange(1, 100, 0.01)
# A =  function(x) #k0(x)
# B = np.zeros_like(A)
# t0 = time.time()
# for i in range(len(x)):
#     B[i] = k0_inv(A[i], function, functionP)
#     # print(f"{x[i]}: {A[i]} -> {B[i]}")
# t1 = time.time()
# print(f"{(t1 - t0)} s for {len(A)} => {1000 * (t1 - t0) / len(A)}ms")
# plt.plot(B, x)
# plt.show()
# exit()

def diffuse_decay_step(C, D, lam, dx, dy, dt, bc="neumann"):
    if bc == "neumann":
        Cpad = np.pad(C, 1, mode="edge")
    elif bc == "dirichlet":
        Cpad = np.pad(C, 1, mode="constant", constant_values=0.0)
    else:
        raise ValueError("bc must be 'neumann' or 'dirichlet'")

    lap = (
        (Cpad[2:, 1:-1] - 2*Cpad[1:-1, 1:-1] + Cpad[:-2, 1:-1]) / dx**2
        + (Cpad[1:-1, 2:] - 2*Cpad[1:-1, 1:-1] + Cpad[1:-1, :-2]) / dy**2
    )

    # allow scalar or spatially varying D
    C_new = C + dt * ( (D * lap) - lam * C )
    return C_new


# --- example usage ---
nx, ny = 128, 128
dx = dy = 1.0
D = 1.0
lam = 0.1
dt = 0.01  # choose small enough for stability
t = 100
C_fixed = np.zeros((ny, nx))
C_inject = np.zeros_like(C_fixed)
q = 10.0
C_inject[ny // 2, nx // 2] = q

for _ in range(int(t/dt)):
    C_fixed[ny // 2, nx // 2] = q  # point mass
    C_fixed = diffuse_decay_step(C_fixed, D, lam, dx, dy, dt, bc="neumann")

    C_inject[ny // 2, nx // 2] += q * dt / (dx*dy)
    # C_inject[ny // 2, int(nx * 2 / 3)] += q * dt / (dx*dy)
    C_inject = diffuse_decay_step(C_inject, D, lam, dx, dy, dt, bc="neumann")

# sigma2 = np.exp(0.5) * (1 - lam)
# C2_fixed = np.zeros_like(C)

xx = np.arange(0, nx, dtype=np.float32)
yy = np.arange(0, ny + 1, dtype=np.float32)
X, Y = np.meshgrid(xx, yy)
D2 = (X - nx/2) ** 2 + (Y - ny / 2) ** 2
r = np.sqrt(D2)

# C2 = np.exp(-5 * np.sqrt(D2) / (D / lam))
C2_fixed =  q * k0(np.sqrt(D2) * np.sqrt(lam/D)) / k0(0.1 * np.sqrt(lam/D)) # (1/(4 * np.pi * D * t)) * np.exp(-D2 / (4 * D * t) - (lam * t))
# C2_inject = q/(2 * np.pi * D) * k0(np.sqrt(D2) * np.sqrt(lam/D))
inte = 0
for s in np.arange(dt, t, dt):
    inte += (1/s) * np.exp(-D2/(4*D*s) - lam * s) * dt
C2_inject = (q/(4*D*np.pi)) * inte

values = C_inject[ny // 2, nx // 2:]

results = np.zeros_like(values)
s = np.arange(dt, t, dt)  # outside the i-loop
for i in range(len(values)):
    def f(r):
        kernel = np.exp(-(r * r) / (4 * D * s) - lam * s) / s
        return (q / (4 * np.pi * D)) * np.sum(kernel) * dt - values[i]

    def fp(r):
        kernel = np.exp(-(r * r) / (4 * D * s) - lam * s) / (s ** 2)
        return -(q * r / (8 * np.pi * D * D)) * np.sum(kernel) * dt
    r0 = i
    try:
        r  = newton(f, r0, fprime=fp, tol=0.1, maxiter=50)
        print(f"{i}: {r}")
        results[i] = r
    except:
        print(f"Not found for {i}")
plt.plot(results)
plt.show()
# exit()

# plt.imshow(np.abs(C2 - C))
# C now contains the field after 50 steps
plt.imshow(C_inject, vmin=0, vmax=1)

# plt.plot(C[ny // 2, nx // 2:])
# plt.plot(C2_fixed[ny // 2, nx // 2:])
# plt.yscale("log")
# plt.plot(np.sqrt(np.maximum(0, -2 * sigma2**2 * np.log(C[ny // 2, nx // 2:] * sigma2 * np.sqrt(2 * np.pi)))))
# plt.plot(np.maximum(0, -2 * sigma2**2 * np.log(C[ny // 2, nx // 2:] * sigma2 * np.sqrt(2 * np.pi))))

# plt.plot(np.abs(C_fixed[ny // 2, nx // 2:] - C2_fixed[ny // 2, nx // 2:]), label="diff")
# plt.plot(np.abs(C_inject[ny // 2, nx // 2:] - C2_inject[ny // 2, nx // 2:]), label="diff")
# plt.plot(C_fixed[ny // 2, nx // 2:], label="Simulation (fixed)")
# plt.plot(C_inject[ny // 2, nx // 2:], label="Simulation (inject)")
# plt.plot(C2_fixed[ny // 2, nx // 2:], label="Analytic (fixed)")
# plt.plot(C2_inject[ny // 2, nx // 2:], label="Analytic (inject)")
# plt.legend()
plt.show()
"""
dt = 0.001  # choose small enough for stability
C2 = np.zeros((ny, nx))

for _ in range(int(t/dt)):
    C2[ny // 2, nx // 2] = 1.0  # point mass
    C2 = diffuse_decay_step(C2, D, lam, dx, dy, dt, bc="neumann")
# C now contains the field after 50 steps
print(np.max(np.abs(C - C2)))
plt.imshow(np.abs(C - C2))
plt.show()
"""
exit(0)




def gaussian(A:np.ndarray, sigma:float):
    H, W = A.shape
    cx, cy = W//2, H//2
    r_px = H
    r = int(max(1, round(r_px)))
    x0 = int(max(0, math.floor(cx - r)))
    x1 = int(min(W-1, math.ceil(cx + r)))
    y0 = int(max(0, math.floor(cy - r)))
    y1 = int(min(H-1, math.ceil(cy + r)))
    if x0 >= x1 or y0 >= y1:
        return
    # sigma = 0.5 * float(r_px)
    # if sigma <= 0.0:
    #     return

    xx = np.arange(x0, x1 + 1, dtype=np.float32)
    yy = np.arange(y0, y1 + 1, dtype=np.float32)
    X, Y = np.meshgrid(xx, yy)
    D2 = (X - cx) ** 2 + (Y - cy) ** 2

    # Truncate to radius r_px
    mask = D2 <= (r_px ** 2)
    weights = np.zeros_like(D2, dtype=np.float32)
    # weights[mask] = np.exp(-0.5 * D2[mask] / (sigma * sigma))
    weights[mask] = (1.0/(sigma * np.sqrt(2))) * np.exp(-0.5 * D2[mask] / (sigma * sigma))
    return weights



















# Ellipsoid scaling demos in matplotlib (pyplot)
# Scenario: a=b=c=1, rotation 45° about the X-axis, and target volume Q=1.0
#
# We show three constructions:
# 1) Amplitude-scaled cap: f = λ g, with λ = 3Q / (2π a b c)
# 2) Axis-scaled cap: use c' = 3Q / (2π a b) (so the half-ellipsoid volume equals Q)
# 3) Cap by clipping plane z0' chosen so the ellipsoidal cap volume equals Q
#
# Notes:
# - We evaluate the rotated ellipsoid in global coordinates via the quadratic form Q = R A R^T.
# - For (1), we scale *height along the local z' normal*; geometrically this moves each cap point
#   r to r' = r + (λ - 1) (n⋅r) n, where n is the unit normal of the clipping plane in global coords.
# - For (2), we change c (in the local frame), then rotate and clip at z' >= 0.
# - For (3), we keep a=b=c=1, solve for h so the volume equals Q, set the plane to z' = c - h,
#   then rotate and keep z' >= z0'.
#
# Implementation choices for plotting:
# - We sample a grid in XY that bounds the projected footprint; we mask outside points.
# - Each method is plotted in its own figure (no subplots), per instructions.
# - We avoid setting explicit colors or styles.

import numpy as np
import matplotlib.pyplot as plt

# ---------- Utility functions ----------

def rotation_x(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[1.0, 0.0, 0.0],
                     [0.0, c,   -s ],
                     [0.0, s,    c ]])

def build_Q(a, b, c, R):
    A = np.diag([1.0/a**2, 1.0/b**2, 1.0/c**2])
    return R @ A @ R.T

def footprint_M(Q):
    Q11, Q22, Q33 = Q[0,0], Q[1,1], Q[2,2]
    Q12, Q13, Q23 = Q[0,1], Q[0,2], Q[1,2]
    v = np.array([[Q13], [Q23]])
    M = Q[:2,:2] - (1.0/Q33) * (v @ v.T)
    return M

def z_upper_from_Q(Q, x, y, translate_z = 0.0):
    Q11, Q22, Q33 = Q[0,0], Q[1,1], Q[2,2]
    Q12, Q13, Q23 = Q[0,1], Q[0,2], Q[1,2]
    A = Q33
    B = 2.0*(Q13*x + Q23*y)
    C = Q11*x**2 + 2.0*Q12*x*y + Q22*y**2 - 1.0
    disc = B*B - 4.0*A*C
    if A == 0 or disc < 0:
        return np.nan
    return (-B + np.sqrt(disc)) / (2.0*A) - translate_z

def ellipse_semi_axes_from_M(M):
    # For ellipse { [x y] M [x y]^T <= 1 }, semi-axes are 1/sqrt(eigs)
    w, _ = np.linalg.eigh(M)
    # Guard: numerical safety in case of tiny negatives
    w = np.maximum(w, 1e-12)
    return 1.0/np.sqrt(w)  # lengths along principal axes

def grid_over_ellipse(M, n=200, pad=1.05):
    s1, s2 = ellipse_semi_axes_from_M(M)
    smax = float(max(s1, s2)) * pad
    xs = np.linspace(-smax, smax, n)
    ys = np.linspace(-smax, smax, n)
    X, Y = np.meshgrid(xs, ys)
    # Mask: inside ellipse -> [x y] M [x y]^T <= 1
    XY = np.stack([X, Y], axis=-1)
    quad = XY[...,0]* (M[0,0]*XY[...,0] + M[0,1]*XY[...,1]) + \
           XY[...,1]* (M[1,0]*XY[...,0] + M[1,1]*XY[...,1])
    inside = quad <= 1.0 + 1e-9
    return X, Y, inside

def sample_cap(Q, z0_prime, R, n=200):
    """
    Build the (x,y,z) arrays for the upper shell of the rotated ellipsoid,
    clipped by the plane z' >= z0' (in the local frame), and translate
    the capped ellipsoid so it *touches* the clipping plane.
    """
    M = footprint_M(Q)
    X, Y, inside = grid_over_ellipse(M, n=n)
    Z = np.full_like(X, np.nan, dtype=float)

    # Plane normal in global coordinates is the third column of R
    nvec = R[:,2]  # unit vector

    # Evaluate z (upper root) where inside the footprint
    idx = np.where(inside)
    for i, j in zip(idx[0], idx[1]):
        x, y = X[i,j], Y[i,j]
        z = z_upper_from_Q(Q, x, y, z0_prime)
        if not np.isnan(z):
            # Check clipping plane: n^T r >= z0'
            r = np.array([x, y, z])
            if nvec @ r >= 0:
                Z[i, j] = z
            # if nvec @ r >= z0_prime - 1e-9:
            #     Z[i,j] = z

    return X, Y, Z, nvec


def ellipsoidal_cap_height_for_volume(Q_target, a, b, c, tol=1e-10, maxit=100):
    """
    Solve for cap height h (in local z') such that the volume of the ellipsoidal cap equals Q_target.
    Cap height h is measured from the clipping plane to the top: 0 <= h <= c.
    Volume formula: V(h) = (π a b / 3) * (h^2 / c) * (3c - h).
    """
    def V(h):
        return (np.pi * a * b / 3.0) * (h*h / c) * (3.0*c - h)

    # Check range
    Vmax = (2.0/3.0) * np.pi * a * b * c  # half-ellipsoid
    if Q_target < 0 or Q_target > Vmax + 1e-12:
        raise ValueError("Q is out of range for a single cap: must be in [0, (2/3)πabc].")

    # Bisection on [0, c]
    lo, hi = 0.0, c
    for _ in range(maxit):
        mid = 0.5*(lo + hi)
        vmid = V(mid)
        if abs(vmid - Q_target) < tol:
            return mid
        if vmid < Q_target:
            lo = mid
        else:
            hi = mid
    return 0.5*(lo+hi)

# ---------- Parameters ----------

a = b = c = 1.0
Q_target = 0.1
angle = 45.0
theta = np.deg2rad(angle)
R = rotation_x(theta)
Qmat = build_Q(a, b, c, R)


# ---------- Method 1: amplitude-scaled cap (f = λ g) ----------

lambda_amp = 3.0 * Q_target / (2.0 * np.pi * a * b * c)
X1, Y1, Z1, nvec = sample_cap(Qmat, z0_prime=0.0, R=R, n=220)
Z_plane = (nvec[0]*X1 - nvec[1]*Y1) / nvec[2]

# Move points along the local-normal direction by (λ - 1) times their signed height above the plane.
# Signed height above the plane equals n·r because the plane passes through the origin and n is unit.
if np.any(~np.isnan(Z1)):
    S = nvec[0]*X1 + nvec[1]*Y1 + nvec[2]*Z1  # signed height above plane
    X1_scaled = X1 + (lambda_amp - 1.0) * S * nvec[0]
    Y1_scaled = Y1 + (lambda_amp - 1.0) * S * nvec[1]
    Z1_scaled = Z1 + (lambda_amp - 1.0) * S * nvec[2]
else:
    X1_scaled = X1.copy()
    Y1_scaled = Y1.copy()
    Z1_scaled = Z1.copy()
# --- Apply the max() so the surface never goes below the clipping plane ---
Z1_clipped = Z_plane.copy()
Z1_clipped[Z1_scaled == Z1_scaled] = np.maximum(Z1_scaled, Z_plane)[Z1_scaled == Z1_scaled]

# ---------- Method 2: axis-scaled cap (set c' so half-volume equals Q) ----------

c2 = 3.0 * Q_target / (2.0 * np.pi * a * b)
Qmat2 = build_Q(a, b, c2, R)
X2, Y2, Z2, _ = sample_cap(Qmat2, z0_prime=0.0, R=R, n=220)
# --- Apply the max() so the surface never goes below the clipping plane ---
Z2_clipped = Z_plane.copy()
Z2_clipped[Z2 == Z2] = np.maximum(Z2, Z_plane)[Z2 == Z2]

# ---------- Method 3: cap by clipping plane z0' (volume of partial ellipsoid equals Q) ----------

h = ellipsoidal_cap_height_for_volume(Q_target, a, b, c)
z0_prime = c - h  # plane is z' = z0' (keep z' >= z0')
X3, Y3, Z3, _ = sample_cap(Qmat, z0_prime=z0_prime, R=R, n=220)

# --- Apply the max() so the surface never goes below the clipping plane ---
Z3_clipped = Z_plane.copy()
Z3_clipped[Z3 == Z3] = np.maximum(Z3, Z_plane)[Z3 == Z3]


# ---------- Print derived quantities ----------

print("Derived parameters:")
print(f"lambda (amplitude scaling) = {lambda_amp:.6f}")
print(f"c' (axis scaling for half-cap volume Q) = {c2:.6f}")
print(f"h (cap height for volume Q) = {h:.6f}")
print(f"z0' (clipping plane in local z') = {z0_prime:.6f}")
print(f"n (plane normal in global coords) = {nvec}")

# ---------- Plotting (one figure per method; no subplots, no explicit colors) ----------

x_lims = [-1.2, 1.2]
y_lims = [-1.2, 1.2]
z_lims = [-1, 1]

# Method 1 figure
fig1 = plt.figure()
ax1 = fig1.add_subplot(131, projection='3d')
ax1.plot_surface(X1, Y1, Z1_clipped, linewidth=0, antialiased=True)
ax1.set_title("Method 1:\namplitude-scaled cap")
ax1.set_xlim(x_lims)
ax1.set_ylim(y_lims)
ax1.set_zlim(z_lims)
ax1.view_init(elev=0, azim=0, roll=0)

fig2 = fig1
fig3 = fig1
# Method 2 figure
# fig2 = plt.figure()
ax2 = fig2.add_subplot(132, projection='3d')
ax2.plot_surface(X2, Y2, Z2_clipped, linewidth=0, antialiased=True)
ax2.set_title("Method 2:\naxis-scaled cap")
ax2.set_xlim(x_lims)
ax2.set_ylim(y_lims)
ax2.set_zlim(z_lims)
ax2.view_init(elev=0, azim=0, roll=0)

# Method 3 figure
# fig3 = plt.figure()
ax3 = fig3.add_subplot(133, projection='3d')
ax3.plot_surface(X3, Y3, Z3_clipped, linewidth=0, antialiased=True)
# ax3.plot_surface(X3, Y3, Z3 - .1, linewidth=0, antialiased=True)
ax3.set_title("Method 3:\ncap by clipping plane $z_0$'")
ax3.set_xlim(x_lims)
ax3.set_ylim(y_lims)
ax3.set_zlim(z_lims)
ax3.view_init(elev=0, azim=0, roll=0)

# plt.suptitle(f"Q={Q_target}, rotation={angle}°, a = b = c = 1.0")
plt.show()








































# import math
#
# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.animation as animation
#
#
# def rotation_matrix(theta, phi, degrees=True):
#     """
#     Construct rotation matrix that rotates the ellipsoid's z-axis to the (theta, phi) direction.
#     theta: azimuthal angle (rotation around z-axis)
#     phi: inclination from z-axis (tilt)
#     """
#     if degrees:
#         theta = np.radians(theta)
#         phi = np.radians(phi)
#
#     # Rotation around z (azimuth)
#     Rz = np.array([
#         [np.cos(theta), -np.sin(theta), 0],
#         [np.sin(theta), np.cos(theta), 0],
#         [0, 0, 1]
#     ])
#
#     # Rotation around y (inclination)
#     Ry = np.array([
#         [np.cos(phi), 0, np.sin(phi)],
#         [0, 1, 0],
#         [-np.sin(phi), 0, np.cos(phi)]
#     ])
#
#     # Final rotation: R = Rz * Ry
#     return Rz @ Ry
#
#
# def ellipsoid_coefficients(a, b, c, R):
#     """
#     Returns the coefficients of the quadratic form of the rotated ellipsoid.
#     The ellipsoid is defined by x.T @ Q @ x = 1, where Q = R @ A @ R.T
#     """
#     A = np.diag([1 / a ** 2, 1 / b ** 2, 1 / c ** 2])
#     Q = R @ A @ R.T  # The rotated quadric matrix
#     return Q
#
# def ellipsoid_roots(x, y, Q, R, z0=0.0, inverse_cut=False):
#     """
#     Given (x, y), return (z_lower_clipped, z_upper_clipped), where clipping occurs
#     at the plane z' = z0 in the local frame, intersecting a global vertical ray.
#     """
#     Qxx, Qxy, Qxz = Q[0, 0], Q[0, 1], Q[0, 2]
#     Qyy, Qyz = Q[1, 1], Q[1, 2]
#     Qzz = Q[2, 2]
#
#     # Solve quadratic in z for vertical line at (x, y)
#     A = Qzz
#     B = 2 * (Qxz * x + Qyz * y)
#     C = Qxx * x**2 + 2 * Qxy * x * y + Qyy * y**2 - 1
#
#     D = B**2 - 4 * A * C
#     if D < 0:
#         return (np.nan, np.nan)
#
#     sqrt_D = np.sqrt(D)
#     z1 = (-B + sqrt_D) / (2 * A)
#     z2 = (-B - sqrt_D) / (2 * A)
#
#     z_upper = max(z1, z2)
#     z_lower = min(z1, z2)
#
#
#     # Compute intersection of vertical line with the cutting plane z' = z0
#     # Find z such that (R^T @ [x, y, z])[2] == z0
#     # => row 3 of R^T dotted with [x, y, z] == z0
#     # => R[2,0] * x + R[2,1] * y + R[2,2] * z_cut = z0
#     # Solve for z_cut
#     v = R.T[2]  # third row of R^T == third column of R # Don't know why I have to translate there...
#     if abs(v[2]) < 1e-5:
#         # Avoid division by near-zero: plane is vertical
#         return (z_lower, z_upper)  # no intersection with vertical ray
#
#     z_cut = (z0 - v[0]*x - v[1]*y) / v[2]
#
#     if inverse_cut:
#         z_upper, z_lower, z_cut = -z_lower, -z_upper, -z_cut
#
#     # Decide whether to replace upper or lower based on orientation
#     if v[2] > 0:
#         if z_cut > z_upper:
#             z_lower, z_upper = np.nan, np.nan
#         else:
#             # local z' points upward in global space — keep upper, clip lower
#             z_lower = max(z_cut, z_lower)
#             # local z' points downward in global space — keep lower, clip upper
#             z_upper = max(z_cut, z_upper)
#     elif z_cut < z_upper:
#         if z_cut < z_lower:
#             z_lower, z_upper = np.nan, np.nan
#         else:
#             # local z' points upward in global space — keep upper, clip lower
#             z_lower = min(z_cut, z_lower)
#             # local z' points downward in global space — keep lower, clip upper
#             z_upper = min(z_cut, z_upper)
#
#     if inverse_cut:
#         z_upper, z_lower, z_cut = -z_lower, -z_upper, -z_cut
#
#     return (z_lower, z_upper)
#
# def normalize(v):
#     norm=np.linalg.norm(v)
#     if norm==0:
#         norm=np.finfo(v.dtype).eps
#     return v/norm
#
# def rotation_from_target(target_x: np.array, target_y: np.array, target_z: np.array):
#     return np.column_stack([target_x, target_y, target_z])
#
#
# def h_plus(x, y, Q, R, z0, inverse_cut=False):
#     return ellipsoid_roots(x, y, Q, R, z0, inverse_cut)[1]
# def h_minus(x, y, Q, R, z0, inverse_cut=False):
#     return ellipsoid_roots(x, y, Q, R, z0, inverse_cut)[0]
# def h(x, y, Q, R, z0, inverse_cut=False):
#     plus = h_plus(x, y, Q, R, z0, inverse_cut)
#     minus = h_minus(x, y, Q, R, z0, inverse_cut)
#     return 0 if plus != plus else plus - minus
#
#
# def height(x, y):
#     octaves = 5
#     lacun = 2.0
#     gain = 1.5
#     return (sum([(math.cos(x * lacun**i * 0.1) + math.sin(15 + y * lacun**i)) / (gain**i) for i in range(octaves)]) / (2.0 * (1 - (1/gain)**(octaves + 1)) * 2.0)) * 0.5 + 0.5
#
# # Ellipsoid and rotation setup
# a, b, c = 0.250, 0.250, 0.150  # semi-axes
# theta = 0               # fixed azimuthal rotation
#
# # Grid
# x_range = np.linspace(-1.1, 1.1, 100)
# y_range = np.linspace(-1.1, 1.1, 100)
# X, Y = np.meshgrid(x_range, y_range)
#
#
# # --- Animation Setup ---
# fig = plt.figure(figsize=(10, 7))
# ax = fig.add_subplot(111, projection='3d')
# surf1 = [None]
# surf2 = [None]
# surf3 = [None]
#
# cx, cy = 0, 0
#
# speedX = 0*0.05
# speedY = 0*0.02
#
# shift1 = 20
# shift2 = -25
#
# def update(frame):
#     global surf1, surf2, surf3, cx, cy
#     ax.clear()
#     z0 = -0.0
#     phi = frame  # degrees
#
#     cx += speedX
#     cy += speedY
#
#     def rotation_and_ellipsoid_at(cx, cy):
#         center_height = height(cx, cy)
#         dx = (center_height - height(cx + 0.01, cy)) / 0.01
#         dy = (center_height - height(cx, cy + 0.01)) / 0.01
#         target_X = normalize(np.array([1, 0, -dx])) # - dx / np.linalg.norm(dx)
#         target_Y = normalize(np.array([0, 1, -dy])) # - dy / np.linalg.norm(dy)
#         target_Z = np.cross(target_X, target_Y)
#
#         R = rotation_from_target(target_X, target_Y, target_Z)
#         # R = rotation_matrix(theta, phi)
#         Q = ellipsoid_coefficients(a, b, c, R)
#         return R, Q
#
#     center_height = height(cx, cy)
#     R1, Q1 = rotation_and_ellipsoid_at(cx + shift1, cy)
#     R2, Q2 = rotation_and_ellipsoid_at(cx + shift2, cy)
#     R, Q = rotation_and_ellipsoid_at(cx, cy)
#     # Z1 = np.array([[h_minus(x, y, Q, R, z0, False) for x, y in zip(x_row, y_row)] for x_row, y_row in zip(X, Y)])
#     # Z2 = np.array([[h_plus(x, y, Q, R, z0, True) for x, y in zip(x_row, y_row)] for x_row, y_row in zip(X, Y)])
#     Z3 = np.array([[h(x, y, Q1, R1, z0, False) for x, y in zip(x_row, y_row)] for x_row, y_row in zip(X, Y)])
#     Z4 = np.array([[h(x, y, Q2, R2, z0, True) for x, y in zip(x_row, y_row)] for x_row, y_row in zip(X, Y)])
#
#     # surf1[0] = ax.plot_surface(X, Y, Z1, cmap='viridis', edgecolor='none', alpha=0.9)
#     # surf2[0] = ax.plot_surface(X, Y, Z2, cmap='plasma', edgecolor='none', alpha=0.9)
#
#     # ax.plot_surface(X, Y, Z3)
#
#
#     Surface = np.array([[height(cx + x * 2.0, cy + y * 2.0) for x, y in zip(x_row, y_row)] for x_row, y_row in zip(X, Y)])
#     # Z1[Z1 != Z1] = 0
#     # Z2[Z2 != Z2] = 0
#     Z3[Z3 != Z3] = -np.inf
#     Z4[Z4 != Z4] = -np.inf
#
#     # surf3[0] = ax.plot_surface(X, Y, np.maximum(Surface, (Z2 - Z1)) - 1.5)
#     surf3[0] = ax.plot_surface(X, Y, (Surface - np.roll(Z3, shift1, axis=0) + np.roll(Z4, shift2, axis=0)) - 1.5)
#
#     # surf3[0] = ax.plot_surface(X, Y, np.maximum(Surface, np.roll(Z3, shift1, axis=0) + height(cx + shift1 * 2.0, cy)) - 1.5)
#     # -- Cutting Plane: z' = z0 in local frame, transformed to global
#     plane_size = 1.2
#     x_plane = np.linspace(-plane_size, plane_size, 20)
#     y_plane = np.linspace(-plane_size, plane_size, 20)
#     X_plane, Y_plane = np.meshgrid(x_plane, y_plane)
#     Z_plane_local = np.full_like(X_plane, z0)
#
#     # Stack local coordinates (x', y', z0)
#     plane_points_local = np.stack([X_plane, Y_plane, Z_plane_local], axis=-1)
#     plane_points_global = plane_points_local @ R.T  # shape (N, N, 3)
#
#     # Extract transformed (x, y, z)
#     X_cut = plane_points_global[:, :, 0]
#     Y_cut = plane_points_global[:, :, 1]
#     Z_cut = plane_points_global[:, :, 2]
#
#     # Plot the cutting plane
#     # ax.plot_surface(X_cut, Y_cut, Z_cut + center_height - 1.5, color='gray', alpha=0.4, edgecolor='none')
#
#
#     ax.set_title(f"Rotated Ellipsoid — phi = {phi:.1f}°")
#     ax.set_xlabel("x")
#     ax.set_ylabel("y")
#     ax.set_zlabel("z")
#     ax.set_xlim(-1.5, 1.5)
#     ax.set_ylim(-1.5, 1.5)
#     ax.set_zlim(-1.5, 1.5)
#
# ani = animation.FuncAnimation(fig, update, frames=np.linspace(0, 360, 60), interval=10)
# # update(45)
# plt.show()