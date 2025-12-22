# Poisson-disk + single Delaunay (pre-removal) + image-weighted border indexing
# + height assignment + 512x512 render.
# Requires: numpy, scipy, pillow (optional: matplotlib)
# Default number of points
nbPoints = 100000

from dataclasses import dataclass, field
import math
import random
from typing import List, Set, Tuple, Optional

import numpy as np
from scipy.spatial import Delaunay

try:
    from PIL import Image
    _HAS_PIL = True
except Exception:
    _HAS_PIL = False

try:
    import matplotlib.pyplot as plt
    _HAS_MPL = True
except Exception:
    _HAS_MPL = False


@dataclass
class Point:
    x: float
    y: float
    index: int = 0           # will be set by indexing
    height: float = 0.0      # will be set as index / max_index
    weight: float = 1          # cached pixel weight in [1..255]; 0 means discarded
    neighbors: Set[int] = field(default_factory=set)  # indices in the ORIGINAL full list

    @property
    def pos(self) -> Tuple[float, float]:
        return (self.x, self.y)


# -----------------------
# Poisson-disk generation
# -----------------------

def bridson_poisson_disk(width: float, height: float, r: float, k: int = 30,
                         seed: Optional[int] = None) -> List[Tuple[float, float]]:
    if seed is not None:
        random.seed(seed)

    cell_size = r / math.sqrt(2)
    grid_w = int(math.ceil(width / cell_size))
    grid_h = int(math.ceil(height / cell_size))
    grid = [None] * (grid_w * grid_h)

    def grid_index(px: float, py: float) -> Tuple[int, int]:
        return int(px // cell_size), int(py // cell_size)

    def fits(px: float, py: float) -> bool:
        if not (0 <= px < width and 0 <= py < height):
            return False
        gi, gj = grid_index(px, py)
        for j in range(max(0, gj - 2), min(grid_h, gj + 3)):
            for i in range(max(0, gi - 2), min(grid_w, gi + 3)):
                q = grid[j * grid_w + i]
                if q is None:
                    continue
                qx, qy = q
                dx, dy = qx - px, qy - py
                if dx*dx + dy*dy < r*r:
                    return False
        return True

    p0 = (random.random() * width, random.random() * height)
    samples = [p0]
    active = [p0]
    gi, gj = grid_index(*p0)
    grid[gj * grid_w + gi] = p0

    while active:
        i = random.randrange(len(active))
        base = active[i]
        placed = False
        for _ in range(k):
            rad = random.uniform(r, 2*r)
            ang = random.uniform(0.0, 2.0*math.pi)
            px = base[0] + rad * math.cos(ang)
            py = base[1] + rad * math.sin(ang)
            if fits(px, py):
                samples.append((px, py))
                active.append((px, py))
                gi, gj = grid_index(px, py)
                grid[gj * grid_w + gi] = (px, py)
                placed = True
                break
        if not placed:
            active.pop(i)

    return samples


def generate_poisson_points_target_count(
    n_target: int,
    domain_size: Tuple[float, float] = (1.0, 1.0),
    k: int = 30,
    seed: Optional[int] = 1234,
    max_adjust_iters: int = 6
) -> List[Tuple[float, float]]:
    width, height = domain_size
    area = width * height
    alpha = 0.55
    r = math.sqrt(area / (n_target * math.pi * alpha))

    points = []
    for t in range(max_adjust_iters):
        pts = bridson_poisson_disk(width, height, r, k=k, seed=None if seed is None else seed + t)
        if len(pts) >= n_target:
            points = pts
            break
        r *= 0.85
    else:
        points = pts

    if len(points) > n_target:
        idxs = np.random.default_rng(seed).choice(len(points), size=n_target, replace=False)
        points = [points[i] for i in idxs]
    elif len(points) < n_target:
        print(f"[warn] Generated {len(points)} points (< target {n_target}).")

    return points


# -----------------------
# Delaunay neighbors (on FULL set)
# -----------------------

def build_delaunay_neighbors(points_xy: List[Tuple[float, float]]) -> List[Set[int]]:
    if len(points_xy) < 3:
        return [set() for _ in points_xy]
    pts = np.asarray(points_xy, dtype=float)
    tri = Delaunay(pts)
    neighbors = [set() for _ in range(len(points_xy))]
    for (a, b, c) in tri.simplices:
        neighbors[a].update((b, c))
        neighbors[b].update((a, c))
        neighbors[c].update((a, b))
    return neighbors


# -----------------------
# Image utilities
# -----------------------

def load_grayscale_image(path: str) -> np.ndarray:
    if not _HAS_PIL:
        raise RuntimeError("Pillow (PIL) is required. Install with `pip install pillow`.")
    img = Image.open(path).convert("L")
    return np.array(img, dtype=np.uint8)

def pixel_value_at(img_u8: np.ndarray, x: float, y: float) -> int:
    h, w = img_u8.shape
    ix = int(round(x * (w - 1)))
    iy = int(round(y * (h - 1)))
    ix = max(0, min(w - 1, ix))
    iy = max(0, min(h - 1, iy))
    return int(img_u8[iy, ix])


# --------------------------------------------
# Border selection + indexing using ORIGINAL neighbors
# --------------------------------------------

def assign_indices_with_image(
    points: List[Point],          # neighbors refer to full set (pre-removal)
    weight_img: np.ndarray,       # uint8 (H,W) values 0..255
    seed: Optional[int] = 42,
    ensure_all_indexed: bool = True
):
    """
    - Keep if pixel>0 (cache weight=pixel 1..255), remove if pixel==0.
    - Initial Border: kept nodes with at least one removed neighbor (using ORIGINAL neighbors).
    - Build kept list and **restrict original neighbors** to kept nodes (no new triangulation).
    - Border nodes get indices 0..|B|-1 (sorted by (x,y)).
    - While Border not empty: pick p from Border with probability ~ p.weight;
      pick a random unindexed neighbor q; if found, index q and add to Border;
      else remove p from Border.
    - Optionally reseed for leftover disconnected components (same weighted rule).
    - Finally, assign height = index / max_index for all kept nodes.
    """
    rng = random.Random(seed)

    # Evaluate mask + cache weights
    pix_vals = [pixel_value_at(weight_img, p.x, p.y) for p in points]
    keep_mask = [pv > 0 for pv in pix_vals]
    removed_mask = [not k for k in keep_mask]
    for p, pv, keep in zip(points, pix_vals, keep_mask):
        p.weight = int(pv if keep else 0)

    # Initial border in original index space
    initial_border_old_idx = set()
    for i, (p, keep) in enumerate(zip(points, keep_mask)):
        if not keep:
            continue
        if any(removed_mask[j] for j in p.neighbors):
            initial_border_old_idx.add(i)

    # Build kept list with remapping and RESTRICT original neighbors
    old_to_new = {}
    kept_points: List[Point] = []
    for old_i, (p, keep) in enumerate(zip(points, keep_mask)):
        if keep:
            old_to_new[old_i] = len(kept_points)
            kept_points.append(
                Point(
                    x=p.x, y=p.y,
                    index=-1, height=0.0, weight=(p.weight / 255.0),
                    neighbors=set()  # fill next
                )
            )

    # Restrict neighbors: for each kept node, keep only neighbors that are kept, remap to new ids
    for old_i, keep in enumerate(keep_mask):
        if not keep:
            continue
        new_i = old_to_new[old_i]
        for nb_old in points[old_i].neighbors:
            if keep_mask[nb_old]:
                kept_points[new_i].neighbors.add(old_to_new[nb_old])

    # Translate initial border to new ids
    border_new_idx = sorted(old_to_new[i] for i in initial_border_old_idx)

    # Assign indices to border nodes deterministically
    border_new_idx.sort(key=lambda i: (kept_points[i].x, kept_points[i].y))
    for idx, pi in enumerate(border_new_idx):
        kept_points[pi].index = idx

    current_index = len(border_new_idx)

    # Weighted choice helper
    def weighted_choice(border_list: List[int]) -> int:
        ws = [kept_points[i].weight for i in border_list]
        return rng.choices(border_list, weights=ws, k=1)[0]

    # Growth phase (weighted by p.weight)
    border_list = list(border_new_idx)
    border_set = set(border_new_idx)

    while border_list:
        p_idx = weighted_choice(border_list)
        unmarked = [q for q in kept_points[p_idx].neighbors if kept_points[q].index < 0]
        if unmarked:
            q_idx = rng.choice(unmarked)
            kept_points[q_idx].index = current_index
            current_index += 1
            if q_idx not in border_set:
                border_set.add(q_idx)
                border_list.append(q_idx)
        else:
            border_list.remove(p_idx)
            border_set.discard(p_idx)

    # Optionally finish leftovers (disconnected components)
    if ensure_all_indexed:
        leftovers = [i for i, p in enumerate(kept_points) if p.index < 0]
        while leftovers:
            seed_i = rng.choice(leftovers)
            kept_points[seed_i].index = current_index
            current_index += 1
            border_list = [seed_i]
            border_set = {seed_i}
            while border_list:
                p_idx = weighted_choice(border_list)
                unmarked = [q for q in kept_points[p_idx].neighbors if kept_points[q].index < 0]
                if unmarked:
                    q_idx = rng.choice(unmarked)
                    kept_points[q_idx].index = current_index
                    current_index += 1
                    if q_idx not in border_set:
                        border_set.add(q_idx)
                        border_list.append(q_idx)
                else:
                    border_list.remove(p_idx)
                    border_set.discard(p_idx)
            leftovers = [i for i, p in enumerate(kept_points) if p.index < 0]

    # Heights = index / max_index
    if kept_points:
        max_idx = max(p.index for p in kept_points)
        if max_idx <= 0:
            for p in kept_points:
                p.height = 0.0
        else:
            inv = 1.0 / float(max_idx)
            for p in kept_points:
                p.height = p.index * inv

    removed_count = sum(1 for k in keep_mask if not k)
    return kept_points, removed_count


# -----------------------
# Rendering 512x512
# -----------------------
def render_points_to_image(points: List[Point], out_path: str, size: int = 512, radius: int = 10):
    """
    Render points to a size×size grayscale image. Each point paints a filled
    disk of 'radius' pixels with value = round(height*255). Overlaps keep the max.
    """
    canvas = np.zeros((size, size), dtype=np.uint8)

    # Precompute squared radius once
    r2 = radius * radius

    # We'll create a local 2D mask per point using 2D broadcast (no extra axes)
    for p in points:
        val = int(max(0, min(255, round(p.height * 255.0))))
        cx = int(round(p.x * (size - 1)))
        cy = int(round(p.y * (size - 1)))
        cx = max(0, min(size - 1, cx))
        cy = max(0, min(size - 1, cy))

        # Local bounding box, clamped to image
        y_min = max(0, cy - radius)
        y_max = min(size, cy + radius + 1)
        x_min = max(0, cx - radius)
        x_max = min(size, cx + radius + 1)

        # Build 2D coordinate grids for the subregion
        # ys: (H,1), xs: (1,W) -> broadcast to (H,W)
        ys = np.arange(y_min, y_max)[:, None]
        xs = np.arange(x_min, x_max)[None, :]

        mask = (xs - cx) * (xs - cx) + (ys - cy) * (ys - cy) <= r2

        sub = canvas[y_min:y_max, x_min:x_max]
        # keep the brightest value on overlaps
        sub[mask] = np.maximum(sub[mask], val)
        canvas[y_min:y_max, x_min:x_max] = sub

    if not _HAS_PIL:
        raise RuntimeError("Pillow (PIL) is required to save images. Install with `pip install pillow`.")
    Image.fromarray(canvas, mode="L").save(out_path)


# -----------------------
# Demo / glue
# -----------------------

def main(
    n: int = nbPoints,
    seed: Optional[int] = 1234,
    image_path: Optional[str] = None,   # grayscale mask 0..255
    out_path: str = "indexed_height.png",
    plot: bool = False
):
    # 1) Poisson-disk points
    pts_xy = generate_poisson_points_target_count(n, domain_size=(1.0, 1.0), seed=seed)

    # 2) Delaunay neighbors on the FULL set (pre-removal)
    full_neighbors = build_delaunay_neighbors(pts_xy)

    # 3) Wrap points with original neighbors
    points: List[Point] = []
    for i, (x, y) in enumerate(pts_xy):
        points.append(Point(x=x, y=y, index=0, height=0.0, weight=1, neighbors=set(full_neighbors[i])))

    # 4) Load grayscale mask
    if image_path is None:
        print("[info] No image provided; using an all-white mask (weights=255 everywhere).")
        weight_img = np.full((512, 512), 255, dtype=np.uint8)
    else:
        weight_img = load_grayscale_image(image_path)

    # 5) Indexing (image-based, using ORIGINAL neighbors restricted to kept nodes)
    kept_points, removed_count = assign_indices_with_image(
        points, weight_img=weight_img, seed=2025, ensure_all_indexed=True
    )

    # 6) Render 512x512
    render_points_to_image(kept_points, out_path=out_path, size=512, radius=2)
    print(f"Saved image to: {out_path}")

    # Some stats
    degrees = [len(p.neighbors) for p in kept_points]
    avg_deg = float(np.mean(degrees)) if degrees else 0.0
    print(f"Kept points: {len(kept_points)} (removed {removed_count}) | "
          f"Avg (restricted) degree: {avg_deg:.2f}")
    if kept_points:
        idxs = sorted(p.index for p in kept_points)
        print(f"Index range: [{idxs[0]} .. {idxs[-1]}] (unique={len(set(idxs))})")
        print(f"Height range: [{min(p.height for p in kept_points):.3f} .. {max(p.height for p in kept_points):.3f}]")

    # Optional: visualize triangulation (still legal, but note neighbors shown here
    # are from the kept set after restriction of original neighbors).
    if plot and _HAS_MPL and len(kept_points) >= 3:
        pts = np.array([[p.x, p.y] for p in kept_points])
        tri = Delaunay(pts)  # purely for plotting triangles over kept points
        plt.triplot(pts[:, 0], pts[:, 1], tri.simplices, linewidth=0.6)
        sc = plt.scatter(pts[:, 0], pts[:, 1], s=6, c=[p.height for p in kept_points])
        plt.colorbar(sc, label="height (index / max_index)")
        plt.title("Kept points (colored by height)")
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    return kept_points


if __name__ == "__main__":
    # Example:
    # main(n=nbPoints, seed=1234, image_path="mask.png", out_path="indexed_height.png", plot=False)
    main(n=nbPoints, seed=1234, image_path="heightmap_mask.png", out_path="indexed_height.png", plot=True)
