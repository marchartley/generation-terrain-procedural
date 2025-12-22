import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401, needed for 3D projection
from matplotlib import cm

# ---------- 1. Define the 1D generating function ----------
_randomProfileCurve = [
    0.8627757352941178, 0.857671568627451, 0.8425245098039216, 0.775888480392157,
    0.6893382352941178, 0.6479779411764706, 0.5974264705882353, 0.5606617647058822,
    0.5392156862745099, 0.508578431372549, 0.4986213235294116, 0.48406862745098034,
    0.4810049019607843, 0.4779411764705883, 0.4756433823529411, 0.47411151960784303,
    0.4595588235294117, 0.4564950980392156, 0.44117647058823517, 0.431985294117647,
    0.42585784313725483, 0.315563725490196, 0.2481617647058823, 0.015318627450980338,
    0.015318627450980338, 0.015, 0.01, 0.001, 0.0, 0.0
]

profile_z = np.array(_randomProfileCurve)
profile_r = np.linspace(0.0, 1.0, len(profile_z))

def z_of_r(r):
    """
    Linear interpolation of z as a function of radius r in [0, 1].
    r can be a scalar or a NumPy array.
    """
    return np.interp(r, profile_r, profile_z) * 0.5

# Parameter along radius
r_min, r_max = 0.0, 1.0
n_r = 200
r = np.linspace(r_min, r_max, n_r)

# Vertical limits based on profile
z_min, z_max = profile_z.min(), profile_z.max()
z_color_min = z_of_r(0.0)
z_color_max = z_of_r(1.0)

# ---------- 2. Prepare the figure ----------
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

ax.set_xlim(-r_max, r_max)
ax.set_ylim(-r_max, r_max)
ax.set_zlim(z_min, z_max)
ax.set_box_aspect([1, 1, (z_max - z_min) / r_max])

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

ax.xaxis.pane.set_visible(False)
ax.yaxis.pane.set_visible(False)
ax.zaxis.pane.set_visible(False)
#
# ax.set_xticks([])
# ax.set_yticks([])
# ax.set_zticks([])
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])

# ax.set_title("Surface of Revolution Animation (Angle-based Stretching)")

# Plot the generating curve (profile) at theta = 0, in the X–Z plane (Y = 0)
z_curve = z_of_r(r)
# ax.plot(r, np.zeros_like(r), z_curve, linewidth=2)

# Plot the rotation axis (the vertical axis through the center)
ax.plot([0, 0], [0, 0], [z_min, z_max], linestyle="--", linewidth=1)

# ---------- 3. Animation & deformation ----------
n_theta = 180          # total angular resolution for a full revolution
n_frames = 100         # how many animation frames

# Predefine all theta values for the full 0..2π range
theta_full = np.linspace(0.0, 2 * np.pi, n_theta)

# Stretching center and amplitude
r0 = 0.0   # center in radius where stretching is "anchored"
a = 0.3    # amplitude of stretch (30%)

def deform_radius_by_angle(R_param, theta_subset):
    """
    Deform the radius field R_param around r0 with a factor that depends
    only on the angle theta, not on the frame.

    This means each angular slice gets a fixed deformation once it exists.
    """
    # Map angle θ ∈ [0, 2π] to a local 'time' t ∈ [0, 1]
    t_local = theta_subset / (2 * np.pi)        # shape (k,)

    # Stretch factor per angle (breathing along the revolution)
    # Each angle 'remembers' this factor forever.
    k_theta = 1.0 + a * np.sin(4 * 2 * np.pi * t_local)  # shape (k,)
    # Broadcast to match R_param (n_r, k)
    # k_theta = k_theta[np.newaxis, :]  # shape (1, k)

    # Shift around center, scale, shift back
    # R_shift = R_param - r0
    # R_def = r0 + k_theta * R_shift
    A = R_param
    B = 3*R_param*R_param - 2*R_param*R_param*R_param
    R_def = (A * (k_theta) + B * (1 - k_theta))
    # Clamp to [0, 1]
    R_def = np.clip(R_def, 0.0, 1.0) ** (1 + k_theta)
    return R_def

# Store a global reference to the current surface
current_surface = None

def init():
    return []

def update(frame):
    global current_surface

    # Remove previous surface if it exists
    if current_surface is not None:
        current_surface.remove()

    # Global fraction of revolution completed (0 -> 1)
    frac = frame / (n_frames - 1)

    # Determine how many theta slices are 'revealed' so far
    max_idx = int(frac * (n_theta - 1)) + 1  # at least 1 slice
    theta_subset = theta_full[:max_idx]      # angles actually drawn

    # Parameter grid for current angles
    R_param, Theta = np.meshgrid(r, theta_subset, indexing="ij")  # shapes (n_r, max_idx)

    # Deform radii based only on their angle (frozen once angle is chosen)
    R_def = deform_radius_by_angle(R_param, theta_subset)

    # Geometry
    Z = z_of_r(R_def)
    X = R_param * np.cos(Theta)
    Y = R_param * np.sin(Theta)

    norm = (Z - z_color_min) / (z_color_max - z_color_min)
    norm = np.clip(norm, 0.0, 1.0)
    colors = cm.gray(norm)  # same shape as Z, RGBA per vertex

    # Draw new surface and save reference
    current_surface = ax.plot_surface(X, Y, Z, alpha=0.7, rcount=100, ccount=100, facecolors=colors, linewidth=0, edgecolor='none')

    return []

# ---------- 4. Create and show the animation ----------
ani = animation.FuncAnimation(
    fig,
    update,
    frames=n_frames,
    init_func=init,
    blit=False,
    interval=50,  # ms between frames
)

ani.save("/home/marc/surface_of_revolution_angle_stretch.mp4", writer="ffmpeg", fps=30)  # optional
plt.show()
