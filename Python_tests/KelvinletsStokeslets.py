from functools import lru_cache
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.interpolate import griddata
from scipy.ndimage import map_coordinates
from imageio.v2 import imread
from PIL import Image


# Paramètres globaux
mu = 0.5           # Viscosité dynamique (ou module de cisaillement)
poisson = 0.25
epsilon = 0.25      # Paramètre de régularisation
force = np.array([1.0, 0.0])  # Force appliquée (vers la droite)
x0 = np.array([0.0, 0.0])     # Point d'application de la force

scale_force = 2.0
twist_force = 2.0

# Grille 2D pour visualisation
width = 200
x_range = np.linspace(-2, 2, width)
y_range = np.linspace(-2, 2, width)
X, Y = np.meshgrid(x_range, y_range)
points = np.stack([X.ravel(), Y.ravel()], axis=1)

x_range = np.linspace(-2, 2, 20)
y_range = np.linspace(-2, 2, 20)
quiver_X, quiver_Y = np.meshgrid(x_range, y_range)
quiver_points = np.stack([quiver_X.ravel(), quiver_Y.ravel()], axis=1)


# display_fn = lambda ax, X, Y, points, field: ax.imshow(compute_divergence_field(points, field).reshape(Y.shape), extent=[-2, 2, -2, 2], cmap='seismic')
# display_fn = lambda ax, X, Y, points, field: ax.imshow(abs(compute_curl_field(points, field)).reshape(Y.shape), extent=[-2, 2, -2, 2], cmap='Reds', vmin=0, vmax=2) #, cmap='seismic')
display_fn = lambda ax, X, Y, points, field: ax.imshow(deformed_img(img, field, points, X, Y).reshape(Y.shape), extent=[-2, 2, -2, 2], cmap='gray')

# img = imread("lena.jpg")
# img = np.array(Image.fromarray(img).resize(X.shape[::-1]))

def create_grid_image(height, width, W=2, S=20):
    img = np.ones((height, width), dtype=np.uint8) * 255  # white background

    # Draw vertical bars
    for x in range(S//2, width, S):
        img[:, x:x + W] = 0

    # Draw horizontal bars
    for y in range(S//2, height, S):
        img[y:y + W, :] = 0

    return img

img = create_grid_image(height=width, width=width, W=2, S=30)


def compute_kelvinlet_grab(x, x0, F, mu, epsilon, poisson = 0.3):
    r = x - x0
    r_norm = np.linalg.norm(r)
    r_eps = np.sqrt(r_norm**2 + epsilon**2)
    I = np.eye(2)
    if r_norm < 1e-10:
        return np.zeros(2)
    term1 = ((3 - 4 * poisson) * I) / r_eps
    term2 = np.outer(r, r) / r_eps**3
    u = (term1 + term2) @ F / (8 * np.pi * mu * (1 - poisson))
    return u

def compute_stokeslet_grab(x, x0, F, mu, epsilon):
    r = x - x0
    r_norm_sq = np.sum(r**2)
    denom = (r_norm_sq + epsilon**2)**(3/2)
    I = np.eye(2)
    term1 = (r_norm_sq + 2 * epsilon**2) * I
    term2 = np.outer(r, r)
    u = ((term1 + term2) @ F) / (8 * np.pi * mu * denom)
    return u

# Définition des primitives supplémentaires pour Kelvinlets et Stokeslets

def kelvinlet_twist(x, x0, torque, mu, epsilon):
    """Twist Kelvinlet - rotation autour de x0"""
    r = x - x0
    r_norm_sq = np.sum(r**2)
    r_eps = np.sqrt(r_norm_sq + epsilon**2)
    if np.linalg.norm(r) < 1e-10:
        return np.zeros(2)
    perp = np.array([-r[1], r[0]])  # rotation de r de 90°
    factor = 1 / r_eps**3
    u = torque * perp * factor / (8 * np.pi * mu)
    return u

def kelvinlet_scale(x, x0, strength, mu, epsilon):
    """Scale Kelvinlet - expansion/contrainte autour de x0"""
    r = x - x0
    r_norm_sq = np.sum(r**2)
    r_eps = np.sqrt(r_norm_sq + epsilon**2)
    if np.linalg.norm(r) < 1e-10:
        return np.zeros(2)
    u = strength * r / r_eps**3 / (8 * np.pi * mu)
    return u

def kelvinlet_pinch(x, x0, direction, magnitude, mu, epsilon):
    """Pinch Kelvinlet - compression dans une direction donnée"""
    r = x - x0
    r_norm_sq = np.sum(r**2)
    r_eps = np.sqrt(r_norm_sq + epsilon**2)
    if np.linalg.norm(r) < 1e-10:
        return np.zeros(2)
    proj = np.dot(r, direction)
    u = -magnitude * proj * direction / r_eps**3 / (8 * np.pi * mu)
    return u


def stokeslet_twist(x, x0, torque, mu, epsilon):
    """Vortex de Stokes - champ de rotation induit (analogue à twist Kelvinlet)"""
    r = x - x0
    r_norm_sq = np.sum(r**2)
    if r_norm_sq < 1e-10:
        return np.zeros(2)
    perp = np.array([-r[1], r[0]])
    denom = (r_norm_sq + epsilon**2)**(3/2)
    u = torque * perp / (8 * np.pi * mu * denom)
    return u

def stokeslet_scale(x, x0, strength, mu, epsilon):
    """Source de Stokes - expansion isotrope (analogue à scale Kelvinlet)"""
    r = x - x0
    r_norm_sq = np.sum(r**2)
    if r_norm_sq < 1e-10:
        return np.zeros(2)
    r_eps = (r_norm_sq + epsilon**2)**(3/2)
    u = strength * r / (8 * np.pi * mu * r_eps)
    return u

def stokeslet_pinch(x, x0, direction, magnitude, mu, epsilon):
    """Champ compressif directionnel - analogue à pinch Kelvinlet"""
    r = x - x0
    r_norm_sq = np.sum(r**2)
    if r_norm_sq < 1e-10:
        return np.zeros(2)
    proj = np.dot(r, direction)
    r_eps = (r_norm_sq + epsilon**2)**(3/2)
    u = -magnitude * proj * direction / (8 * np.pi * mu * r_eps)
    return u


# Champs associés
def compute_field(points, field_fn):
    return np.array([field_fn(p) for p in points])

def compute_magnitude_field(points, vectors):
    magnitudes = np.linalg.norm(vectors, axis=1)
    return magnitudes

def compute_divergence_field(points, vector_field):
    # Infer grid dimensions from unique x and y coordinates
    x_unique = np.unique(points[:, 0])
    y_unique = np.unique(points[:, 1])
    nx, ny = len(x_unique), len(y_unique)

    # Reshape the field
    U = vector_field[:, 0].reshape(ny, nx)
    V = vector_field[:, 1].reshape(ny, nx)

    dx = x_unique[1] - x_unique[0]
    dy = y_unique[1] - y_unique[0]

    dU_dx = np.gradient(U, dx, axis=1)
    dV_dy = np.gradient(V, dy, axis=0)

    divergence = dU_dx + dV_dy
    return divergence

def compute_curl_field(points, vector_field):
    x_unique = np.unique(points[:, 0])
    y_unique = np.unique(points[:, 1])
    nx, ny = len(x_unique), len(y_unique)

    U = vector_field[:, 0].reshape(ny, nx)  # u_x
    V = vector_field[:, 1].reshape(ny, nx)  # u_y

    dx = x_unique[1] - x_unique[0]
    dy = y_unique[1] - y_unique[0]

    dV_dx = np.gradient(V, dx, axis=1)
    dU_dy = np.gradient(U, dy, axis=0)

    curl = dV_dx - dU_dy
    return curl

def arrows(field, quiver_points, ax, quiver_X, quiver_Y):
    Uq = griddata(points, field[:, 0], quiver_points, method='linear').reshape(quiver_Y.shape)
    Vq = griddata(points, field[:, 1], quiver_points, method='linear').reshape(quiver_Y.shape)
    # Compute the magnitude of the original (interpolated) vectors
    magnitudes = np.sqrt(Uq ** 2 + Vq ** 2)
    # Flatten arrays for plotting
    Uq_flat = Uq.ravel()
    Vq_flat = Vq.ravel()
    magnitudes_flat = magnitudes.ravel()
    Xq_flat = quiver_X.ravel()
    Yq_flat = quiver_Y.ravel()
    # Normalize magnitudes for coloring
    norm = Normalize(vmin=np.min(magnitudes_flat), vmax=np.max(magnitudes_flat))
    colors = plt.cm.turbo(norm(magnitudes_flat) * 2)  # Or any other colormap
    # Normalize vectors (limit length to consistent scale)
    epsilon = 1e-8  # to avoid division by zero
    max_length = np.max(np.diff(x_range)) * 0.8  # max arrow length ~80% of a cell
    scaling_factors = np.maximum(magnitudes_flat, epsilon)
    Uq_scaled = Uq_flat / scaling_factors * max_length
    Vq_scaled = Vq_flat / scaling_factors * max_length
    # Plot
    ax.quiver(Xq_flat, Yq_flat, Uq_scaled, Vq_scaled,
              color=colors, angles='xy', scale_units='xy', scale=1)

def deformed_img(img, vector_field, points, X, Y):
    ny, nx = Y.shape
    coords_x = X.copy()
    coords_y = Y.copy()

    # Reshape vector field to 2D grids
    U = vector_field[:, 0].reshape(ny, nx)
    V = vector_field[:, 1].reshape(ny, nx)

    # Compute inverse mapping: where each pixel came from
    map_x = coords_x - U
    map_y = coords_y - V

    # Normalize mapping to image pixel coordinates
    h, w = img.shape[:2]
    xi = ((map_x + 2) / 4) * (w - 1)  # [-2, 2] → [0, w-1]
    yi = ((map_y + 2) / 4) * (h - 1)  # [-2, 2] → [0, h-1]

    # Interpolate image at new coordinates
    if img.ndim == 2:  # grayscale
        warped = map_coordinates(img, [yi.ravel(), xi.ravel()], order=1, mode='reflect')
        return warped.reshape(h, w)
    else:  # RGB
        channels = []
        for c in range(img.shape[2]):
            warped_c = map_coordinates(img[..., c], [yi.ravel(), xi.ravel()], order=1, mode='reflect')
            channels.append(warped_c.reshape(h, w))
        return np.stack(channels, axis=-1)



# Définition des fonctions à appliquer
grab_fn  = lambda p: compute_kelvinlet_grab(p, x0, force, mu=mu, epsilon=epsilon, poisson=poisson)
twist_fn = lambda p: kelvinlet_twist(p, x0, torque=twist_force, mu=mu, epsilon=epsilon)
scale_fn = lambda p: kelvinlet_scale(p, x0, strength=scale_force, mu=mu, epsilon=epsilon)
pinch_fn = lambda p: kelvinlet_pinch(p, x0, direction=np.array(force), magnitude=np.linalg.norm(force), mu=mu, epsilon=epsilon)

# Calculs
grab_field  = compute_field(points, grab_fn)
twist_field = compute_field(points, twist_fn)
scale_field = compute_field(points, scale_fn)
pinch_field = compute_field(points, pinch_fn)

# Affichage
fig, axs = plt.subplots(2, 4, figsize=(15, 5))

fields = [(grab_field, "Grab"), (twist_field, "Twist"), (scale_field, "Scale"), (pinch_field, "Pinch")]
for ax, (field, title) in zip(axs[0], fields):
    display_fn(ax, X, Y, points, field)
    # arrows(field, quiver_points, ax, quiver_X, quiver_Y)
    ax.set_ylim([-2, 2])
    ax.set_xlim([-2, 2])
    ax.set_title(f"Kelvinlet - {title}")
    ax.axis('equal')
    ax.grid(False)

# Application sur la grille
grab_fn  = lambda p: compute_stokeslet_grab(p, x0, force, mu=mu, epsilon=epsilon)
stokes_twist_fn = lambda p: stokeslet_twist(p, x0, torque=twist_force, mu=mu, epsilon=epsilon)
stokes_scale_fn = lambda p: stokeslet_scale(p, x0, strength=scale_force, mu=mu, epsilon=epsilon)
stokes_pinch_fn = lambda p: stokeslet_pinch(p, x0, direction=np.array(force), magnitude=np.linalg.norm(force), mu=mu, epsilon=epsilon)

# Calculs des champs
grab_field  = compute_field(points, grab_fn)
stokes_twist_field = compute_field(points, stokes_twist_fn)
stokes_scale_field = compute_field(points, stokes_scale_fn)
stokes_pinch_field = compute_field(points, stokes_pinch_fn)

fields = [(grab_field, "Grab"), (stokes_twist_field, "Twist"), (stokes_scale_field, "Scale"), (stokes_pinch_field, "Pinch")]
for ax, (field, title) in zip(axs[1], fields):
    display_fn(ax, X, Y, points, field)
    # arrows(field, quiver_points, ax, quiver_X, quiver_Y)
    ax.set_ylim([-2, 2])
    ax.set_xlim([-2, 2])
    ax.set_title(f"Stokeslet - {title}")
    ax.axis('equal')
    ax.grid(False)

plt.tight_layout()
plt.show()
