# from __future__ import annotations
import copy
import math
import random
import itertools
from typing import List, Union, Tuple, Optional, Any

import matplotlib.backend_bases
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import scipy.spatial
from heightmap_extractor.constraints_tests import check_feasibility, Vector2D, circle_coords, str_distances
import warnings
warnings.filterwarnings("ignore")


class Node:
    currentIDs: List = []  # static var
    id: int

    connectedNodes: List[Tuple['Node', float]]

    def __init__(self, _id=None):
        self.connectedNodes = []
        if _id is None:
            _id = max(Node.currentIDs) + 1 if Node.currentIDs else 0
        self.id = _id
        Node.currentIDs.append(_id)

    def connect(self, node: 'Node', weight: float = 1.0):
        if len([n for n, w in self.connectedNodes if n.id == node.id]) > 0:
            return
        self.connectedNodes.append((node, weight))

    def disconnect(self, node: Union['Node', int]):
        removedId = node.id if isinstance(node, Node) else node
        self.connectedNodes = [n for n in self.connectedNodes if n[0].id != removedId]

    def degree(self) -> int:
        return len(self.connectedNodes)

    def isConnected(self, node: Union['Node', int]):
        return len([n.id for n, w in self.connectedNodes if (isinstance(node, Node) and n.id == node.id) or n.id == node]) > 0

    def __repr__(self):
        return str(self)

    def __str__(self):
        out = "Node #" + str(self.id) + ". " + str(len(self.connectedNodes)) + " neighbor(s)" + \
              ((":" + ", ".join(["#" + str(n.id) + "(w=" + str(round(w, 2)) + ")" for n, w in self.connectedNodes])) if self.connectedNodes else "")
        return out


class Graph:
    nodes: List[Node]

    def __init__(self, nodes: List[Node] = None):
        if nodes is None:
            self.nodes = []
        else:
            self.nodes = nodes

    def findNode(self, node: Union[Node, int]) -> Node:
        foundId = node.id if isinstance(node, Node) else node
        returned = [n for n in self.nodes if n.id == foundId]
        return returned[0] if returned else None

    def findIndex(self, node: Node) -> int:
        returnedIndices = [i for i, n in enumerate(self.nodes) if n.id == node.id]
        return returnedIndices[0] if returnedIndices else None

    def addNode(self, node: Node = None):
        if node is None:
            node = Node()
        self.nodes.append(node)

    def removeNode(self, node: Union[Node, int]):
        removedId = None
        if isinstance(node, Node):
            removedId = node.id
        else:
            removedId = node
        self.nodes = [n for n in self.nodes if n.id != removedId]
        for i in range(len(self.nodes)):
            self.nodes[i].id = i

    def addNeighbors(self, nodeA: Union[Node, int], nodeB: Union[Node, int], weightA: float = 1.0, weightB: float = None):
        weightB = weightA if weightB is None else weightA
        self.findNode(nodeA).connect(self.findNode(nodeB), weightA)
        self.findNode(nodeB).connect(self.findNode(nodeA), weightB)

    def __repr__(self):
        return f"Graph containing {len(self.nodes)} node(s)"

    def __str__(self):
        return f"Graph containing {len(self.nodes)} node(s)" + \
               (":\n" + "\n".join([str(n) for n in self.nodes]) if self.nodes else "")

    def adjacencyMatrix(self, ignoreSelfConnection: bool = False) -> np.ndarray:
        adj = self.weightedAdjacencyMatrix(ignoreSelfConnection)
        adj[adj > 0] = 1
        adj[adj < 0] = 0
        return adj

    def weightedAdjacencyMatrix(self, ignoreSelfConnection: bool = False) -> np.ndarray:
        adj = np.zeros((len(self.nodes), len(self.nodes)))
        for i, node in enumerate(self.nodes):
            connectedNodesIndices = [(self.findIndex(n), weight) for n, weight in node.connectedNodes]
            for j, w in connectedNodesIndices:
                if ignoreSelfConnection and i == j:
                    continue
                adj[i, j] = w
        return adj

    def normalizedAdjacencyMatrix(self) -> np.ndarray:
        invD = self.inverseSqrtDegreeMatrix()
        return np.matmul(np.matmul(invD, self.adjacencyMatrix()), invD)  # = D^{-1/2} A D^{-1/2}

    def degreeMatrix(self) -> np.ndarray:
        degree = np.zeros((len(self.nodes), len(self.nodes)))
        for i in range(len(self.nodes)):
            degree[i, i] = self.nodes[i].degree()
        return degree

    def inverseSqrtDegreeMatrix(self) -> np.ndarray:
        invDegree = np.zeros((len(self.nodes), len(self.nodes)))
        for i in range(len(self.nodes)):
            invDegree[i, i] = 1 / math.sqrt(max(self.nodes[i].degree(), 1))
        return invDegree

    def laplacianMatrix(self) -> np.ndarray:
        return self.degreeMatrix() - self.adjacencyMatrix()

    def normalizedLaplacianMatrix(self) -> np.ndarray:
        return np.identity(len(self.nodes)) - self.normalizedAdjacencyMatrix()

    def transitionMatrix(self) -> np.ndarray:
        return np.matmul(np.linalg.inv(self.degreeMatrix()), self.adjacencyMatrix())
        # return np.matmul(np.linalg.inv(self.degreeMatrix()), self.normalizedAdjacencyMatrix())

    def edgesIndices(self) -> List[Tuple[int, int]]:
        result = []
        for i, node in enumerate(self.nodes):
            neighbors = [n for n, w in node.connectedNodes]
            result += [(i, self.findIndex(n)) for n in neighbors]
        return result

    def edgesIndicesWithWeight(self) -> List[Tuple[int, int, float]]:
        result = []
        for i, node in enumerate(self.nodes):
            neighbors = node.connectedNodes
            result += [(i, self.findIndex(n), w) for n, w in neighbors]
        return result

    def makeUniqueComponent(self):
        if not self.nodes:
            return
        # adjacency = self.adjacencyMatrix(ignoreSelfConnection=True)
        mainComponent = set()  # [self.nodes[0].id]) # + [n.id for n, w in self.nodes[0].connectedNodes])
        for i in range(len(self.nodes)):
            _id = self.nodes[i].id
            if _id not in mainComponent:
                randomIndex = i
                while (len(mainComponent) == 0 and randomIndex == i) or (
                        len(mainComponent) > 0 and randomIndex not in mainComponent):
                    randomIndex = random.randint(0, len(self.nodes) - 1)
                mainComponent.add(_id)
                self.addNeighbors(self.nodes[i], self.nodes[randomIndex])
            mainComponent = mainComponent.union(set([n.id for n, w in self.nodes[i].connectedNodes]))
        return

    def getAllFaces(self) -> List[List[Node]]:
        newNodes = []
        paths = []

        for i in range(len(self.nodes)):
            newNodes.append(self.nodes[i])
            paths.append([self.nodes[i]])

        cycleFound = False
        shortestCycles = []
        while not cycleFound and len(newNodes) > 0:
            newNodes.clear()
            newPaths = []
            paths = list(filter(lambda p: p, paths))
            for iPath, path in enumerate(paths):
                lastNode = path[-1]
                newNodes.append(lastNode)
                for neighbor, w in lastNode.connectedNodes:
                    if neighbor not in path:
                        newNodes.append(neighbor)
                        newPaths.append(path + [neighbor])
                    elif len(path) >= 3 and neighbor == path[0]:
                        # cycleFound = True
                        shortestCycles.append(path)
            paths = newPaths

        # Filter in multiple passes to improve performance
        # First, remove any cycle containing a smaller cycle
        filteredCycles: List = []
        for i in inv_range(len(shortestCycles)):
            isValid = True
            for j in range(i):
                if containsSubarray(shortestCycles[i], shortestCycles[j], checkReverse=False, checkSlidingWindow=False):
                    isValid = False
                    break
            if isValid:
                filteredCycles.append(shortestCycles[i])
        shortestCycles = filteredCycles
        filteredCycles = []
        # Then try again by considering the inverse of the cycle (0-1-2-3 == 3-2-1-0)
        # and also a rolling cycle (0-1-2-3 == 2-3-0-1)
        for i in inv_range(len(shortestCycles)):
            isValid = True
            for j in range(i):
                if containsSubarray(shortestCycles[i], shortestCycles[j], checkReverse=True, checkSlidingWindow=True):
                    isValid = False
                    break
            if isValid:
                filteredCycles.append(shortestCycles[i])

        # shortestCycles = filteredCycles
        # filteredCycles = []
        # for i in inv_range(len(shortestCycles)):
        #     isValid = True
        #     if len(shortestCycles[i]) != 3:
        #         for nodeA in shortestCycles[i]:
        #             for nodeB in shortestCycles[i]:
        #                 for nodeC in shortestCycles[i]:
        #                     if nodeA.id != nodeB.id and nodeA.id != nodeC.id and nodeB.id != nodeC.id:
        #                         if nodeC in [n for n, w in nodeA.connectedNodes] and nodeC in [n for n, w in nodeB.connectedNodes]:
        #                             isValid = False
        #     if isValid:
        #         filteredCycles.append(shortestCycles[i])

        return filteredCycles  # shortestCycles

    def triangulate(self, randomWeights=None):
        if randomWeights is None:
            randomWeights = [1.0, 1.0]
        faces = self.getAllFaces()
        faces = list(filter(lambda face: len(face) > 3, faces))

        for face in faces:
            while len(face) > 3:
                pivotVertex = random.randint(0, len(face))
                self.addNeighbors(face[pivotVertex - 1], face[(pivotVertex + 1) % len(face)], randomWeights[0] + random.random() * (randomWeights[1] - randomWeights[0]))
                face = face[:pivotVertex] + face[pivotVertex+1:]

    def copy(self):
        return copy.deepcopy(self)

    @classmethod
    def randomPlanarGraph(cls, numberOfPoints: int):
        points = [[random.random(), random.random()] for _ in range(numberOfPoints)]
        voro = scipy.spatial.Voronoi(points)

        graph = cls([Node(i) for i in range(numberOfPoints)])
        for idA, idB in voro.ridge_points:
            graph.addNeighbors(idA, idB, compute_dist(voro.points[idA], voro.points[idB]) * 20.0)
        return graph

    @staticmethod
    def gnp_random_connected_graph(n, p):
        """
        Generates a random undirected graph, similarly to an Erdős-Rényi
        graph, but enforcing that the resulting graph is connected
        """
        edges = itertools.combinations(range(n), 2)
        G = nx.Graph()
        G.add_nodes_from(range(n))
        if p <= 0:
            return G
        # if p >= 1:
        #     return nx.complete_graph(n, create_using=G)
        for _, node_edges in itertools.groupby(edges, key=lambda x: x[0]):
            node_edges = list(node_edges)
            random_edge = random.choice(node_edges)
            G.add_edge(*random_edge)
            if p >= 1:
                random_edges = random.sample(node_edges, int(p - 1))
                for e in random_edges:
                    G.add_edge(*e)
            else:
                for e in node_edges:
                    if random.random() < p:
                        G.add_edge(*e)
        return G


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


def inv_range(start: int, end: int = None, step: int = 1):
    if end is None:
        end = start
        start = 0

    current = end
    while current > start:
        current -= step
        yield current


def containsSubarray(originalArray: List, subarray: List, checkReverse: bool = False, checkSlidingWindow: bool = False) -> bool:
    notContained = len([x for x in range(len(originalArray)) if originalArray[x:x + len(subarray)] == subarray]) == 0

    if checkSlidingWindow:
        tmp = originalArray
        for _ in range(len(originalArray)):
            tmp = tmp[1:] + [tmp[0]]
            notContained = notContained and len([x for x in range(len(tmp)) if tmp[x:x + len(subarray)] == subarray]) == 0

    if checkReverse:
        notContained = notContained and not containsSubarray(originalArray,
                                                             list(reversed(subarray)),
                                                             False,
                                                             checkSlidingWindow)
    return not notContained
