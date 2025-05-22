# import numpy as np
# import networkx as nx
# import igraph as ig
# import leidenalg as la
# import matplotlib.pyplot as plt
# from scipy.spatial import Delaunay, distance_matrix
# from PIL import Image  # For loading and processing the image
#
# # Step 1: Load image, convert to grayscale, and sample points
# def sample_points_from_image(image_path, num_points, noise=0.02):
#     # Load the image
#     image = Image.open(image_path).convert("L")  # Convert to grayscale
#     image = np.array(image)  # Convert to numpy array
#
#     # Normalize the image to [0, 1] range
#     image = 1 - image / 255.0
#
#     # Get image dimensions
#     height, width = image.shape
#
#     # Flatten the image and normalize to create a probability distribution
#     prob_dist = image.flatten() / np.sum(image)
#
#     # Sample points based on the probability distribution
#     sampled_indices = np.random.choice(len(prob_dist), size=num_points, p=prob_dist)
#     sampled_points = np.array([(idx % width, idx // width) for idx in sampled_indices])
#
#     # Normalize points to [0, 1] range and add noise
#     sampled_points = sampled_points / np.array([width, height])
#     sampled_points += np.random.normal(scale=noise, size=sampled_points.shape)
#
#     return sampled_points
#
# # Step 2: Generate points from an image
# image_path = "test.png"  # Replace with the path to your image
# num_points = 1000  # Number of points to sample
# points = sample_points_from_image(image_path, num_points, noise=0.0)
# D = 0.04
#
# # Step 3: Compute Delaunay triangulation to create edges
# tri = Delaunay(points)
# edges = set()
# edge_weights = {}
#
# dist_matrix = distance_matrix(points, points)
# for i in range(num_points):
#     for j in range(i + 1, num_points):  # Avoid duplicate edges
#         dist = dist_matrix[i, j]
#         if dist < D:
#             edges.add((i, j))
#             edge_weights[(i, j)] = float(1 / (max(dist, 0.01) ** 1))  # Stronger weights for closer nodes
# # # Iterate over simplices (triangles) in the Delaunay triangulation
# # for simplex in tri.simplices:
# #     # Add edges for each triangle
# #     for i in range(3):
# #         for j in range(i + 1, 3):
# #             a, b = simplex[i], simplex[j]
# #             # if a < b:  # Avoid duplicate edges
# #             edges.add((a, b))
# #             dist = np.linalg.norm(points[a] - points[b])
# #             edge_weights[(a, b)] = 1 # float(1 / (dist ** 2))  # Stronger weights for closer nodes
#
# # Step 4: Convert to NetworkX graph
# G_nx = nx.Graph()
# for i, (x, y) in enumerate(points):
#     G_nx.add_node(i, pos=(x, y))
# for (a, b), weight in edge_weights.items():
#     G_nx.add_edge(a, b, weight=weight)
#
# # Convert to iGraph for Leiden
# edges_list = [(a, b, weight) for (a, b), weight in edge_weights.items()]
# G_ig = ig.Graph(n=len(points))  # Initialize graph with all nodes
# G_ig.add_edges([(a, b) for a, b, _ in edges_list])  # Add edges
# G_ig.es['weight'] = [weight for _, _, weight in edges_list]  # Add edge weights
#
# # Step 5: Run Leiden Algorithm with weighted modularity
# partition = la.find_partition(G_ig, la.CPMVertexPartition, weights=G_ig.es['weight'])
#
# # Step 6: Assign colors to communities
# num_communities = len(set(partition.membership))
# cmap = plt.cm.turbo  # Use a colormap with many colors
# community_colors = {c: cmap(i / num_communities) for i, c in enumerate(set(partition.membership))}
# node_colors = [community_colors[partition.membership[v]] for v in range(len(points))]
#
# # Step 7: Visualization
# plt.figure(figsize=(8, 8))
# pos = {i: (points[i][0], points[i][1]) for i in range(len(points))}
# nx.draw(G_nx, pos, with_labels=False, node_color=node_colors, edge_color='gray', node_size=50)
# plt.title("Leiden Algorithm with Delaunay Triangulation")
# print(partition)
# plt.show()


import open3d as o3d
import numpy as np
from scipy.spatial import cKDTree
import igraph as ig
import leidenalg
import matplotlib.pyplot as plt

# Load the Stanford Bunny model
mesh = o3d.io.read_triangle_mesh("bunny.ply")
mesh.compute_vertex_normals()
points = np.asarray(mesh.vertices)
vertices = np.asarray(mesh.vertices)
normals = np.asarray(mesh.vertex_normals)
# Get triangle faces (indices of vertices)
triangles = np.asarray(mesh.triangles)

# Use a set to store unique undirected edges
edge_set = set()
for tri in triangles:
    i, j, k = tri
    edge_set.update([(min(i, j), max(i, j)),
                     (min(j, k), max(j, k)),
                     (min(k, i), max(k, i))])

edges = list(edge_set)
weights = []

# Use dot product of normals to compute edge weights
for i, j in edges:
    ni, nj = normals[i], normals[j]
    similarity = abs(np.dot(ni, nj))  # -1 to 1
    weight = 1/np.linalg.norm(ni - nj) * (similarity + 1e-3)**1  # normalize to [0,1], add epsilon
    weights.append(weight)

# Create igraph graph
g = ig.Graph(edges=edges, directed=False)
g.es['weight'] = weights

# Run Leiden algorithm
partition = leidenalg.find_partition(g, leidenalg.RBERVertexPartition, weights=g.es['weight'])
print(partition)
# Assign colors to vertices by cluster
num_clusters = len(partition)
colors = plt.get_cmap("tab20")(np.linspace(0, 1, num_clusters))[:, :3]  # RGB only

vertex_colors = np.zeros((len(points), 3))
for cluster_idx, cluster in enumerate(partition):
    for vertex_idx in cluster:
        vertex_colors[vertex_idx] = colors[cluster_idx]

# Apply colors to mesh and visualize
mesh.vertex_colors = o3d.utility.Vector3dVector(vertex_colors)
o3d.visualization.draw_geometries([mesh], window_name='Stanford Bunny - Leiden Clusters')
