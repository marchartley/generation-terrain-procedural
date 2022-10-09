import random
from typing import List, Optional, Tuple, Union

from numpy import ndarray
from scipy.spatial import Voronoi

from Graph import Graph, Node
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib

# matplotlib.use('Qt5Agg')

import scipy.spatial
import numpy as np
import math
from matplotlib.tri.triangulation import Triangulation

from heightmap_extractor.constraints_tests import Vector2D, check_feasibility
import networkx as nx


class GeometricGraph(Graph):
    def __init__(self, nodes: List[Node] = None):
        super().__init__(nodes)
        self.positions = np.empty((len(self.nodes), len(self.nodes)))
        self.positions[:] = np.nan

    @classmethod
    def fromGraph(cls, normalGraph: Graph) -> 'GeometricGraph':
        graph: GeometricGraph = GeometricGraph(normalGraph.nodes)
        return graph

    def addNode(self, node: Node = None):
        super().addNode(node)
        self.positions.resize((len(self.nodes), len(self.nodes)))

    def removeNode(self, node: Union[Node, int]):
        super().removeNode(node)
        self.positions.resize((len(self.nodes), len(self.nodes)))

    def meanErrorOnPosition(self):
        error = 0
        percentError = 0
        mse = 0
        total = 0
        nbConnections = 0
        for iNode, node in enumerate(self.nodes):
            for iNeighbor, (neighbor, weight) in enumerate(node.connectedNodes):
                realDistance = compute_dist(self.positions[self.findIndex(node), :2],
                                            self.positions[self.findIndex(neighbor), :2])
                if weight == 0 or realDistance == 0:
                    continue
                nbConnections += 1
                total += weight
                error += abs(weight - realDistance)
                mse += (weight - realDistance) ** 2
                percentError += max(weight / realDistance, realDistance / weight)
                pass
        return [error / nbConnections, (percentError / nbConnections) - 1, mse / nbConnections]

    def displayPositioningError(self) -> str:
        MAE, MPE, MSE = self.meanErrorOnPosition()
        return f"Result :\nMAE = {round(MAE, 3)}\nMSE = {round(MSE, 3)}\n(~{round(100 * MPE, 2)}% error per node)\n" \
               f"{len(self.getAllEdgeIntersections())} intersections"

    def computeVoronoi(self) -> Optional[Tuple[List[List[float]], ndarray, List[List[np.ndarray]], Voronoi]]:
        used_points = [pos[:2] for i, pos in enumerate(self.positions) if not np.any(np.isnan(pos))]
        if len(used_points) > 2:
            voro = scipy.spatial.Voronoi(used_points)
            vertices_indices, regions_coord = voronoi_finite_polygons_2d(voro, 100)
            areas = []
            for indices in vertices_indices:
                areas.append([regions_coord[idx] for idx in indices])
            return vertices_indices, regions_coord, areas, voro
        else:
            return None

    def computeDelaunay(self, limitedNodes: List[Union[Node, int]] = None) -> Optional[Tuple[List[List[int]],
                                                List[List[int]],
                                                List[float],
                                                List[float],
                                                scipy.spatial.Delaunay]]:
        _limitedNodes = []
        if limitedNodes is None:
            _limitedNodes = [n.id for n in self.nodes]
        else:
            for i in range(len(limitedNodes)):
                if isinstance(limitedNodes[i], Node):
                    _limitedNodes.append(limitedNodes[i].id)
                else:
                    _limitedNodes.append(limitedNodes[i])
        used_points = [pos[:2] for i, pos in enumerate(self.positions) if not np.any(np.isnan(pos))]
        vertices_index = [i for i, pos in enumerate(self.positions) if not np.any(np.isnan(pos))]
        if len(used_points) > 2:
            delaunay = scipy.spatial.Delaunay(used_points)
            triangles = [[vertices_index[index] for index in triangle] for triangle in delaunay.simplices if
                         triangle[0] in _limitedNodes and triangle[1] in _limitedNodes and triangle[2] in _limitedNodes]
            neighbors = [[]] * len(self.positions)
            for i in range(len(delaunay.vertex_neighbor_vertices[0]) - 1):
                if vertices_index[i] not in _limitedNodes:
                    continue
                neighbors[vertices_index[i]] = [vertices_index[n] for n in delaunay.vertex_neighbor_vertices[1][
                                                                           delaunay.vertex_neighbor_vertices[0][i]:
                                                                           delaunay.vertex_neighbor_vertices[0][i + 1]] if vertices_index[n] in _limitedNodes]

            x, y = delaunay.points.T
            edges = np.array([e for e in Triangulation(x, y, delaunay.simplices).edges if e[0] in _limitedNodes and e[1] in _limitedNodes])
            tri_lines_x = np.insert(x[edges], 2, np.nan, axis=1).ravel()
            tri_lines_y = np.insert(y[edges], 2, np.nan, axis=1).ravel()
            return neighbors, triangles, tri_lines_x.tolist(), tri_lines_y.tolist(), delaunay
        else:
            return None

    def getSpatialConnectionsBetweenNodes(self) -> Tuple[int, int, int]:
        wrongConnections = 0
        missingConnections = 0
        correctConnections = 0
        delaunayResults = self.computeDelaunay()
        if delaunayResults:
            delaunayNeighbors, _, _, _, _ = delaunayResults
            for i, node in enumerate(self.nodes):
                correct = 0
                missing = 0
                for realNeighbor, w in node.connectedNodes:
                    if realNeighbor.id in delaunayNeighbors[i]:
                        correct += 1
                    else:
                        missing += 1
                wrongConnections += len(delaunayNeighbors[i]) - correct
                correctConnections += correct
                missingConnections += missing

        missingConnections //= 2
        wrongConnections //= 2
        correctConnections //= 2
        return correctConnections, wrongConnections, missingConnections

    def displayGraphInfo(self) -> str:
        MAE, MPE, MSE = self.meanErrorOnPosition()
        nbNodes = len(self.positions)
        nbInstances = len([pos for pos in self.positions if not np.any(np.isnan(pos))])

        correctConnections, wrongConnections, missingConnections = self.getSpatialConnectionsBetweenNodes()
        output = f"MAE  = {round(MAE, 2)}\nMSE  = {round(MSE, 2)}\n%err = {round(100 * MPE, 2)}%\n" \
                 f"{nbInstances}/{nbNodes} nodes placed\n" \
                 f"{correctConnections} good adjacencies\n" \
                 f"{missingConnections} missing adjacencies\n" \
                 f"{wrongConnections} wrong adjacencies\n"
        return output

    def getAllEdgeIntersections(self) -> List[Tuple[Vector2D, Node, Node, Node, Node]]:
        intersections = []
        for i, nodeA in enumerate(self.nodes):
            for j, (nodeB, w1) in enumerate(nodeA.connectedNodes):
                if nodeA.id < nodeB.id:  # Avoid counting edges twice
                    for ii, otherA in enumerate(self.nodes):
                        if i > ii:
                            for jj, (otherB, w2) in enumerate(otherA.connectedNodes):
                                if otherA.id < otherB.id:  # Avoid counting edges twice
                                    inter = line_intersection(self.positions[nodeA.id], self.positions[nodeB.id],
                                                              self.positions[otherA.id], self.positions[otherB.id],
                                                              0.001)
                                    if inter:
                                        intersections.append(
                                            (Vector2D(inter[0], inter[1]), nodeA, nodeB, otherA, otherB))
        return intersections

    mousePressEventCID = None
    mouseMoveEventCID = None
    mouseReleaseEventCID = None
    movingPoint = None
    GUI_Text: matplotlib.pyplot.Text = None

    def draw(self,
             withLabels: bool = False,
             withWeights: bool = False,
             withVoronoi: bool = False,
             withDelaunay: bool = False,
             withCirclesAround: bool = False,
             selectedNodes: List[Node] = None,
             title: str = ""):

        Xs = [pos[0] for pos in self.positions]
        Ys = [pos[1] for pos in self.positions]
        if selectedNodes is None:
            selectedNodes = self.nodes
        for i in range(len(selectedNodes)):
            if isinstance(selectedNodes[i], int):
                selectedNodes[i] = self.findNode(selectedNodes[i])

        if withCirclesAround:
            weights = [min([w for n, w in node.connectedNodes]) / 2 if node.connectedNodes else np.nan for node in
                       self.nodes if node in selectedNodes]
            nbVertexOnCircle = 20
            for i, weight in enumerate(weights):
                plt.plot([Xs[i] + math.cos(a * 2 * math.pi / nbVertexOnCircle) * weight for a in
                          range(nbVertexOnCircle + 1)],
                         [Ys[i] + math.sin(a * 2 * math.pi / nbVertexOnCircle) * weight for a in
                          range(nbVertexOnCircle + 1)])

        if self.mousePressEventCID is not None:
            plt.gcf().canvas.mpl_disconnect(self.mousePressEventCID)
            plt.gcf().canvas.mpl_disconnect(self.mouseMoveEventCID)
            plt.gcf().canvas.mpl_disconnect(self.mouseReleaseEventCID)
            self.mousePressEventCID = self.mouseMoveEventCID = self.mouseReleaseEventCID = None

        if self.GUI_Text is not None:
            self.GUI_Text.remove()

        if withVoronoi:
            voronoiResults = self.computeVoronoi()
            if voronoiResults:
                _, _, _, voro = voronoiResults
                scipy.spatial.voronoi_plot_2d(voro, show_points=False, show_vertices=False, ax=plt.gca())

        plt.scatter(Xs, Ys, c="r")
        for i, node in enumerate(self.nodes):
            if node not in selectedNodes:
                continue
            for _j, neighbor in enumerate(node.connectedNodes):
                j = self.findIndex(neighbor[0])
                if neighbor[0] in selectedNodes:
                    plt.plot([Xs[i], Xs[j]], [Ys[i], Ys[j]], linewidth=2.0, color="green")

        if withDelaunay:
            delaunayResults = self.computeDelaunay(limitedNodes = selectedNodes)
            if delaunayResults:
                _, _, tri_lines_x, tri_lines_y, delaunay = delaunayResults
                plt.plot(tri_lines_x, tri_lines_y, linewidth=1.0, color="red", label="Geometric adjacency")

                def onClick(event):
                    mouse = [event.xdata, event.ydata]
                    closest = math.inf
                    for _i, pos in enumerate(self.positions):
                        if not np.any(np.isnan(pos)):
                            dist = compute_dist(pos[:2], mouse)
                            if dist < closest:
                                closest = dist
                                self.movingPoint = _i

                def onMove(event):
                    if self.movingPoint is not None:
                        mouse = [event.xdata, event.ydata]
                        self.positions[self.movingPoint, :2] = mouse
                        plt.clf()
                        self.draw(withLabels=withLabels, withWeights=withWeights, withVoronoi=withVoronoi,
                                  withDelaunay=withDelaunay, title="Manual editing")
                        plt.draw()

                def onRelease(_event):
                    self.movingPoint = None

                if self.mouseMoveEventCID is None:
                    self.mousePressEventCID = plt.gcf().canvas.mpl_connect('button_press_event', onClick)
                    self.mouseMoveEventCID = plt.gcf().canvas.mpl_connect('motion_notify_event', onMove)
                    self.mouseReleaseEventCID = plt.gcf().canvas.mpl_connect('button_release_event', onRelease)

        if withLabels:
            for i in range(len(self.nodes)):
                if self.nodes[i] in selectedNodes:
                    plt.annotate(str(self.nodes[i].id), (Xs[i], Ys[i]))
        if withWeights:
            for i, node in enumerate(self.nodes):
                if node in selectedNodes:
                    for _j, (neighbor, weight) in enumerate(node.connectedNodes):
                        j = self.findIndex(neighbor)
                        if j not in selectedNodes:
                            continue
                        difference = weight - compute_dist([Xs[i], Ys[i]], [Xs[j], Ys[j]])
                        sDiff = ("+" if difference > 0 else "") + str(round(difference, ndigits=2))
                        plt.annotate(str(round(weight, ndigits=2)) + "\n(" + sDiff + ")",
                                     ((Xs[i] + Xs[j]) / 2, (Ys[i] + Ys[j]) / 2))

        if withDelaunay:
            line1, = plt.plot([], [], linewidth=1.0, color="red", label="Geometric adjacency")
            line2, = plt.plot([], [], color="green", linewidth=2.0, label="Graph adjacency")
            plt.legend(handles=[line1, line2])

        plt.gca().set_ylabel("Y")
        plt.gca().set_xlabel("X")
        plt.gca().set_aspect('equal')
        plt.title(title)

        plt.subplots_adjust(left=.1, right=.7)
        self.GUI_Text = plt.figtext(.75, 0.5, self.displayGraphInfo(), wrap=True)
        plt.gcf()

    def isPerfectlyPositioned(self):
        correctConnections, wrongConnections, missingConnections = self.getSpatialConnectionsBetweenNodes()
        return wrongConnections == 0 and missingConnections == 0


def compute_dist(a, b=None) -> float:
    distance = 0
    if b is None:
        b = [0] * len(a)
    if len(a) != len(b):
        raise Exception("The two vectors have different sizes")
    nbDims = min(len(a), len(b))
    for i in range(nbDims):
        distance += (a[i] - b[i]) ** 2
    return math.sqrt(distance)


def voronoi_finite_polygons_2d(vor: scipy.spatial.Voronoi, radius: float = None) -> Tuple[
    List[List[float]], np.ndarray]:
    """
    From https://gist.github.com/Sklavit/e05f0b61cb12ac781c93442fbea4fb55
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.
    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.
    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.
    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max() * 2

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1]  # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:, 1] - c[1], vs[:, 0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)


def line_intersection(P11, P12, P21, P22, epsilon: float = 0.0) -> Optional[List[float]]:
    if isinstance(P11, Vector2D):
        P11 = P11.asarray()
    if isinstance(P12, Vector2D):
        P12 = P12.asarray()
    if isinstance(P21, Vector2D):
        P21 = P21.asarray()
    if isinstance(P22, Vector2D):
        P22 = P22.asarray()
    t = ((P11[0] - P21[0]) * (P21[1] - P22[1]) - (P11[1] - P21[1]) * (P21[0] - P22[0])) / (
            (P11[0] - P12[0]) * (P21[1] - P22[1]) - (P11[1] - P12[1]) * (P21[0] - P22[0]))
    u = ((P11[0] - P21[0]) * (P11[1] - P12[1]) - (P11[1] - P21[1]) * (P11[0] - P12[0])) / (
            (P11[0] - P12[0]) * (P21[1] - P22[1]) - (P11[1] - P12[1]) * (P21[0] - P22[0]))

    # check if line actually intersect
    if (epsilon != 0.0 and 0 + epsilon < t < 1 - epsilon and 0 + epsilon < u < 1 - epsilon) or (
            epsilon == 0.0 and 0 < t < 1 and 0 < u < 1):
        return [P11[0] + t * (P12[0] - P11[0]), P11[1] + t * (P12[1] - P11[1])]
    else:
        return None


def iterativeGeometryFromLaplacianEigenVectors(graph: Graph, initialPositions: np.ndarray = None) -> np.ndarray:
    if initialPositions is None:
        nodesPositions = np.linalg.eig(graph.normalizedLaplacianMatrix())[1]
    else:
        nodesPositions = initialPositions.copy()
    transition = (np.identity(len(graph.nodes)) + graph.transitionMatrix()) * .5
    nodesPositions = np.matmul(transition, nodesPositions).real
    return nodesPositions


def geometryFromInstancingNodes(graph: Graph, distanceTolerance: float = 1.5, verbose: bool = False) -> np.ndarray:
    nbNodes = len(graph.nodes)
    positions = iterativeGeometryFromInstancingNodes(graph, distanceTolerance=distanceTolerance, verbose=verbose)
    for _ in range(nbNodes):
        positions = iterativeGeometryFromInstancingNodes(graph, positions, distanceTolerance, verbose)
    return positions


def iterativeGeometryFromInstancingNodes(graph: Graph,
                                         initialPositions: np.ndarray = None,
                                         distanceTolerance: float = 1.5,
                                         verbose: bool = True) -> np.ndarray:
    nbNodes = len(graph.nodes)
    if initialPositions is None:
        nodesPositions = np.empty((nbNodes, nbNodes))  # np.zeros((nbNodes, nbNodes))
        nodesPositions[:] = np.nan
    else:
        nodesPositions = initialPositions.copy()

    setNodes = []
    # noNodeSet = True
    for i in range(nbNodes):
        if np.any(np.isnan(nodesPositions[i, :])):
            # allNodesSet = False
            pass
        else:
            # noNodeSet = False
            setNodes.append(i)

    # All nodes are already set
    if len(setNodes) == nbNodes:
        return nodesPositions

    # Nothing is set, just place the first one randomly
    elif len(setNodes) == 0:
        nodesPositions[0] = [0] * nbNodes

    else:
        nodesIndexToInstantiate = -1
        minEntropy = math.inf
        for i in range(nbNodes):
            if i in setNodes:
                continue
            node = graph.nodes[i]
            instantiateNeighbors = [n[0].id for n in node.connectedNodes if n[0].id in setNodes]
            entropy = 1 / max(.1, len(instantiateNeighbors))  # node.degree() / max(1, len(instantiateNeighbors))
            if entropy < minEntropy:
                minEntropy = entropy
                nodesIndexToInstantiate = i
        nodeToInstantiate = graph.findNode(nodesIndexToInstantiate)
        # print(f"Selected node #{nodeToInstantiate.id} with entropy of {minEntropy}")
        possibleAABBox = [math.inf, -math.inf, math.inf, -math.inf]
        for i, (neighbor, w) in enumerate(nodeToInstantiate.connectedNodes):
            if neighbor.id in setNodes:
                possibleAABBox = [
                    min(nodesPositions[neighbor.id, 0] - distanceTolerance * w, possibleAABBox[0]),
                    max(nodesPositions[neighbor.id, 0] + distanceTolerance * w, possibleAABBox[1]),
                    min(nodesPositions[neighbor.id, 1] - distanceTolerance * w, possibleAABBox[2]),
                    max(nodesPositions[neighbor.id, 1] + distanceTolerance * w, possibleAABBox[3]),
                ]
        if possibleAABBox[0] == math.inf:
            possibleAABBox = [-1000, 2000, -1000, 2000]
        else:
            possibleAABBox = [possibleAABBox[0],
                              possibleAABBox[1] - possibleAABBox[0],
                              possibleAABBox[2],
                              possibleAABBox[3] - possibleAABBox[2]]
        numberOfTries = 1000
        positionValidated = False
        newPosition = []
        bestPositionFound = []
        nbCrossingsOnBestPosition = math.inf
        while numberOfTries > 0 and not positionValidated:
            newPosition = [possibleAABBox[0] + random.random() * possibleAABBox[1],
                           possibleAABBox[2] + random.random() * possibleAABBox[3]]
            positionValidated = True
            for i, (neighbor, w) in enumerate(nodeToInstantiate.connectedNodes):
                if neighbor.id in setNodes:
                    neighborPos = nodesPositions[neighbor.id, :2]
                    sqr_distance = (neighborPos[0] - newPosition[0]) ** 2 + (neighborPos[1] - newPosition[1]) ** 2
                    if sqr_distance < (w * w) or \
                            (numberOfTries > 100 and sqr_distance > distanceTolerance * distanceTolerance * w * w) or \
                            (numberOfTries > 100 and sqr_distance < (w * w) / (distanceTolerance * distanceTolerance)):
                        positionValidated = False
                        break

            if positionValidated:
                hypotheticalEdges = []
                for i, (neighbor, w) in enumerate(nodeToInstantiate.connectedNodes):
                    if neighbor.id in setNodes:
                        hypotheticalEdges.append([newPosition] + [nodesPositions[neighbor.id, :2]])

                nbCrossings = 0
                for i, node in enumerate(graph.nodes):
                    if node.id in setNodes:
                        for j, (neighbor, w) in enumerate(node.connectedNodes):
                            if neighbor.id in setNodes:
                                for P1, P2 in hypotheticalEdges:
                                    if line_intersection(P1, P2, nodesPositions[node.id, :2],
                                                         nodesPositions[neighbor.id, :2]):
                                        nbCrossings += 1
                if nbCrossings < nbCrossingsOnBestPosition:
                    nbCrossingsOnBestPosition = nbCrossings
                    bestPositionFound = newPosition
                if nbCrossings > 0:
                    positionValidated = False
            numberOfTries -= 1
        if positionValidated:
            nodesPositions[nodesIndexToInstantiate] = newPosition + [0] * (nbNodes - 2)
            if verbose:
                print(f"{nodeToInstantiate} is set at ({round(newPosition[0], 2)}, {round(newPosition[1], 2)})")
        else:
            nodesPositions[nodesIndexToInstantiate] = bestPositionFound + [0] * (nbNodes - 2)
            if verbose:
                print(f"We had to force {nodeToInstantiate}, "
                      f"creating {nbCrossingsOnBestPosition} crossings in the plane...")
            # if verbose:
            #     print(f"We didn't find a spot for {nodeToInstantiate}...")
    return nodesPositions


def repositionNodesWithForces(graph: GeometricGraph, verbose: bool = False) -> np.ndarray:
    distances = []
    for i, node in enumerate(graph.nodes):
        node_distances = list([None] * len(graph.nodes))
        node_distances[i] = 0.0
        for neighbor, distance in node.connectedNodes:
            node_distances[neighbor.id] = distance
        distances.append(node_distances)

    if verbose:
        if not check_feasibility(distances, verbose=False):
            print("No config possible to get the constraints done :")
            check_feasibility(distances, verbose=True)
        else:
            print("Theoretically we can solve the positioning")

    epsilon = 1e-4
    tries = 0
    maxIterations = 300
    deltaMove = .5
    deltaMoveDamping = 1.

    returningPositions = graph.positions.copy()
    while True:

        newPositions = iterativeRepositionNodesWithForce(graph, deltaMove)

        moves = newPositions - returningPositions
        tries += 1
        if np.all(np.array([compute_dist(move) for move in moves]) < epsilon) or tries >= maxIterations:
            break

        deltaMove *= deltaMoveDamping
        if deltaMove < 1e-4:
            break

        returningPositions = newPositions

    return returningPositions


def iterativeRepositionNodesWithForce(graph: GeometricGraph, deltaMove: float) -> np.ndarray:
    nb_nodes = len(graph.nodes)
    initialPositions = graph.positions.copy()
    nodes = [Vector2D(initialPositions[i, 0], initialPositions[i, 1]) for i in range(nb_nodes)]

    distances = []
    for i, node in enumerate(graph.nodes):
        node_distances = list([None] * nb_nodes)
        node_distances[i] = 0.0
        for neighbor, distance in node.connectedNodes:
            node_distances[neighbor.id] = distance
        distances.append(node_distances)

    for i in range(len(distances)):
        for j in range(len(distances)):
            if distances[i][j] is None:
                distances[i][j] = -1

    for i in range(len(distances)):
        for j in range(len(distances)):
            distances[i][j] = max(distances[i][j], distances[j][i])

    dists = np.zeros((len(nodes), len(nodes)))
    nodes_pairs = []
    for i in range(len(nodes)):
        nodes_pairs.append([])
        for j in range(len(nodes)):
            nodes_pairs[i].append([])
    moves = []
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            nodes_pairs[i][j] = nodes[j] - nodes[i]
            dists[i, j] = (nodes[i] - nodes[j]).norm()

    for i in range(len(nodes)):
        move = Vector2D(0, 0)
        divisor = 0
        moveVectors = []
        for j in range(len(nodes)):
            if distances[i][j] >= 0:
                move += nodes_pairs[i][j].normalized() * (nodes_pairs[i][j].norm() - distances[i][j])
                if (nodes_pairs[i][j].norm() - distances[i][j]) != 0:
                    moveVectors.append([nodes_pairs[i][j].normalized(), (nodes_pairs[i][j].norm() - distances[i][j])])
                divisor += 1
        if divisor > 0:
            move /= divisor
        move *= deltaMove

        moves.append(move)
        nodes[i] += move

    returningPositions = np.zeros_like(initialPositions)
    for i, node in enumerate(nodes):
        returningPositions[i, :2] = [node.x, node.y]
    return returningPositions


def main():
    Node.currentIDs.clear()
    nbNodes = 10
    nbEdgesPerNode = 2
    # nbDimensions = 2
    customDrawing = True

    # nodes = []
    # for i in range(nbNodes):
    #     nodes.append(Node())
    #
    # graph = Graph(nodes)
    # for i in range(nbNodes):
    #     # for j in range(i + 1, nbNodes):
    #     #     graph.addNeighbors(graph.nodes[i], graph.nodes[j])
    #     # if i > 0:
    #     graph.addNeighbors(graph.nodes[i], graph.nodes[(i + 1) % len(graph.nodes)])
    #     # for iEdge in range(nbEdgesPerNode):
    #     #     randomlyChosenNode = i
    #     #     while randomlyChosenNode == i:
    #     #         randomlyChosenNode = random.randint(0, nbNodes - 1)
    #     #     graph.addNeighbors(graph.nodes[i], graph.nodes[randomlyChosenNode], 10 + random.random() * 10)
    #
    # graph.makeUniqueComponent()
    graph = GeometricGraph.randomPlanarGraph(numberOfPoints=nbNodes)
    # graph.triangulate(randomWeights = [10, 20])

    if customDrawing:
        bestTolerance = 1.5
        bestError = math.inf
        for t in [2.0, 1.8, 1.5, 1.3, 1.2, 1.15, 1.05]:
            testingGraph = graph.copy()
            testingGraph.positions = geometryFromInstancingNodes(testingGraph, t)
            errors = len(testingGraph.getAllEdgeIntersections())
            if errors <= bestError:
                bestError = errors
                bestTolerance = t
        distanceTolerance = bestTolerance

        # randomPositions = iterativeGeometryFromLaplacianEigenVectors(graph)
        graph.positions = iterativeGeometryFromInstancingNodes(graph, distanceTolerance=distanceTolerance)

        # angles = [np.vdot(randomPositions[i, :], randomPositions[(i + 1) % nbNodes]) for i in range(nbNodes)]

        plt.ion()
        for iteration in range(nbNodes):
            plt.pause(0.01)
            # randomPositions = iterativeGeometryFromLaplacianEigenVectors(graph, randomPositions)
            graph.positions = iterativeGeometryFromInstancingNodes(graph, graph.positions)
            plt.clf()
            plt.draw()
            graph.draw(withLabels=True, withWeights=False,
                       withVoronoi=True, withDelaunay=True,
                       withCirclesAround=True,
                       title=f"Instancing points ({iteration + 1}/{nbNodes})")

        for iteration in range(100):
            plt.pause(0.05)
            plt.clf()
            plt.cla()
            graph.positions = iterativeRepositionNodesWithForce(graph, .99)
            graph.draw(withLabels=True, withWeights=False,
                       withVoronoi=True, withDelaunay=True,
                       withCirclesAround=True,
                       title=f"Apply forces ({round(100 * (iteration + 1) / 100, 2)}%)")
            plt.draw()
        plt.ioff()

        print(f"Tolerance = {distanceTolerance}\n" + graph.displayPositioningError())
        print(graph.displayPositioningError())
        plt.pause(0.01)
        plt.show(block=True)

    else:
        G = GeometricGraph.gnp_random_connected_graph(nbNodes, nbEdgesPerNode)
        nx.draw_spectral(G, with_labels=True)
        plt.show()


if __name__ == "__main__":
    random.seed(42)
    for _ in range(40):
        main()
