# Side-by-side reaction–diffusion comparison: Gray–Scott, Schnakenberg, Brusselator, Gierer–Meinhardt
# This script simulates each model on a 2D grid and plots the "activator/product" field for each one.
# Notes:
# - Uses explicit Euler + 5-point Laplacian with periodic boundaries.
# - Parameters are chosen to yield visible patterns in a short run; you can tweak steps or params to change looks.
# - No seaborn, one chart per figure rule is followed by making a single composite figure.

import numpy as np
import matplotlib.pyplot as plt

# ---------- Utilities ----------

def laplacian_periodic(Z):
    # 5-point stencil with periodic BCs
    return (
        -4*Z
        + np.roll(Z, 1, axis=0) + np.roll(Z, -1, axis=0)
        + np.roll(Z, 1, axis=1) + np.roll(Z, -1, axis=1)
    )

def simulate_gray_scott(n=128, steps=3000, Du=0.16, Dv=0.08, F=0.035, k=0.065, dt=1.0):
    # Initialize
    u = np.ones((n, n))
    v = np.zeros((n, n))
    # Seed a square in the center
    r = n//10
    u[n//2 - r:n//2 + r, n//2 - r:n//2 + r] = 0.50
    v[n//2 - r:n//2 + r, n//2 - r:n//2 + r] = 0.25
    # Small noise
    u += 0.01*np.random.rand(n, n)
    v += 0.01*np.random.rand(n, n)

    for _ in range(steps):
        Lu = laplacian_periodic(u)
        Lv = laplacian_periodic(v)
        uvv = u * v * v
        u += (Du * Lu - uvv + F * (1 - u)) * dt
        v += (Dv * Lv + uvv - (F + k) * v) * dt
        # clamp to sane range to avoid runaway due to explicit Euler
        u = np.clip(u, 0, 1.5)
        v = np.clip(v, 0, 1.5)

    return u, v

def simulate_schnakenberg(n=128, steps=4000, Du=0.08, Dv=1.0, a=0.2, b=1.3, dt=0.25):
    # steady state (for reference): u* = a + b, v* = b / (a + b)^2
    u = (a + b) * np.ones((n, n))
    v = (b / (a + b)**2) * np.ones((n, n))
    # noise
    u += 0.02*(np.random.rand(n, n) - 0.5)
    v += 0.02*(np.random.rand(n, n) - 0.5)

    for _ in range(steps):
        Lu = laplacian_periodic(u)
        Lv = laplacian_periodic(v)
        uv2 = u*u*v
        u += (Du*Lu + a - u + uv2) * dt
        v += (Dv*Lv + b - uv2) * dt
        u = np.clip(u, 0, 3.0)
        v = np.clip(v, 0, 3.0)
    return u, v

def simulate_brusselator(n=128, steps=3000, Du=0.1, Dv=2.0, A=1.0, B=3.0, dt=0.2):
    # steady state: u* = A, v* = B/A
    u = A * np.ones((n, n))
    v = (B/A) * np.ones((n, n))
    # noise
    u += 0.05*(np.random.rand(n, n) - 0.5)
    v += 0.05*(np.random.rand(n, n) - 0.5)

    for _ in range(steps):
        Lu = laplacian_periodic(u)
        Lv = laplacian_periodic(v)
        uv2 = u*u*v
        u += (Du*Lu + A - (B + 1)*u + uv2) * dt
        v += (Dv*Lv + B*u - uv2) * dt
        u = np.clip(u, 0, 5.0)
        v = np.clip(v, 0, 5.0)
    return u, v

def simulate_gierer_meinhardt(n=128, steps=4000, Da=0.01, Dh=0.4, mu_a=0.02, mu_h=0.06, rho=0.04, dt=0.15):
    # activator a and inhibitor h
    a = 0.5 * np.ones((n, n))
    h = 1.0 * np.ones((n, n))
    # noise
    a += 0.02*(np.random.rand(n, n) - 0.5)
    h += 0.02*(np.random.rand(n, n) - 0.5)

    for _ in range(steps):
        La = laplacian_periodic(a)
        Lh = laplacian_periodic(h)
        # avoid division by zero
        denom = np.maximum(h, 1e-6)
        react_a = (a*a)/denom - mu_a*a + rho
        react_h = a*a - mu_h*h
        a += (Da*La + react_a) * dt
        h += (Dh*Lh + react_h) * dt
        a = np.clip(a, 0, 5.0)
        h = np.clip(h, 1e-6, 10.0)
    return a, h

# ---------- Run all models ----------

np.random.seed(7)  # reproducibility (roughly)

u_gs, v_gs = simulate_gray_scott()
u_sn, v_sn = simulate_schnakenberg()
u_br, v_br = simulate_brusselator()
a_gm, h_gm = simulate_gierer_meinhardt()

# ---------- Plot side-by-side ----------

fig, axs = plt.subplots(2, 2, figsize=(10, 10))

# Show "product/activator" fields that typically display patterns best:
# Gray–Scott: v, Schnakenberg: u or v both interesting (use u), Brusselator: u, Gierer–Meinhardt: a
im0 = axs[0, 0].imshow(v_gs)
axs[0, 0].set_title("Gray–Scott (v)")
axs[0, 0].axis('off')
plt.colorbar(im0, ax=axs[0, 0], fraction=0.046, pad=0.04)

im1 = axs[0, 1].imshow(u_sn)
axs[0, 1].set_title("Schnakenberg (u)")
axs[0, 1].axis('off')
plt.colorbar(im1, ax=axs[0, 1], fraction=0.046, pad=0.04)

im2 = axs[1, 0].imshow(u_br)
axs[1, 0].set_title("Brusselator (u)")
axs[1, 0].axis('off')
plt.colorbar(im2, ax=axs[1, 0], fraction=0.046, pad=0.04)

im3 = axs[1, 1].imshow(a_gm)
axs[1, 1].set_title("Gierer–Meinhardt (activator a)")
axs[1, 1].axis('off')
plt.colorbar(im3, ax=axs[1, 1], fraction=0.046, pad=0.04)

plt.tight_layout()
plt.show()