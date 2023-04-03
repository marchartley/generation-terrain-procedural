import math
import random
import time
from typing import List, Union

from AdjencyMapCreation import getCompatibleNodesForNewFace
from GeometricGraph import geometryFromInstancingNodes, repositionNodesWithForces
from TopoMap import CombinMap, Brin, Graph, GeometricGraph, Node, Vector2D
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import os
import json

from pathfinding.core.grid import Grid as pathfinderGrid
from pathfinding.core.diagonal_movement import DiagonalMovement as pathfinderDiagonalMovement
from pathfinding.finder.bi_a_star import BiAStarFinder as BiAStarFinder

def gridAsRegions(_grid: np.ndarray) -> np.ndarray:
    grid = _grid.astype(int) #.copy()
    for y in range(grid.shape[0]):
        for x in range(grid.shape[1]):
            if grid[x, y] <= 0:
                nodeId = int(abs(grid[x, y]) - 2)
            elif grid[x, y] > 1:
                nodeId = grid[x, y] - 2
            else:
                nodeId = -1
            grid[x, y] = int(nodeId)
    return grid

def displayGrid(grid, nodesPos, sizeX, sizeY):
    for y in range(sizeY):
        for x in range(sizeX):
            nodeId = -1
            for i in range(len(nodesPos)):
                if nodesPos[i] == [x, y]:
                    nodeId = i
                    break
            if grid[x, y] <= 0:
                nodeId = int(abs(grid[x, y]) - 2)
                print(chr(ord("a") + nodeId), end="")
            elif nodeId > -1:
                print(chr(ord("A") + nodeId), end="")
            else:
                print(" ", end="")
        print()

def gridShow(grid: np.ndarray, colors: List[int] = None, fig: plt.Figure = None, tofile: Union[bool, str] = False):
    tmp = gridAsRegions(grid)
    if colors is None:
        colors = [i for i in range(max(np.unique(tmp)) + 1)]
    colored = np.empty(grid.shape, dtype=int)
    for y in range(tmp.shape[0]):
        for x in range(tmp.shape[1]):
            nodeId = tmp[y, -(x + 1)]
            if nodeId >= 0:
                colored[x, y] = colors[nodeId]
            else:
                colored[x, y] = 0

    if tofile:
        plt.imshow(colored, vmin=0)
        plt.savefig(tofile)
    else:
        if fig is None:
            plt.imshow(colored, vmin=0)
        else:
            ax = fig.add_subplot()
            ax.imshow(colored, vmin=0)
        plt.show()

def getRegionsSizes(grid: np.ndarray, maxNodeId: int) -> List[int]:
    # uniques = np.unique(grid)
    sizes: List[int] = []
    for node in range(maxNodeId):
        sizes.append(len(grid[grid == node]))
    return sizes


def dilute(grid: np.ndarray, desiredSizes: List[int]) -> np.ndarray:
    tmp = grid.copy()
    regions = gridAsRegions(grid)
    actualSizes = getRegionsSizes(regions, len(desiredSizes))
    nodesToGrow = np.array(actualSizes) < np.array(desiredSizes)
    for x in range(grid.shape[0]):
        for y in range(grid.shape[1]):
            val = grid[x, y]
            if not nodesToGrow[regions[x, y]]:
                continue
            if val != 1:
                for dx in [-1, 0, 1]:
                    for dy in [-1, 0, 1]:
                        if 0 <= (x + dx) < grid.shape[0] and 0 <= (y + dy) < grid.shape[1] and grid[x + dx, y + dy] == 1:
                            tmp[x + dx, y + dy] = val
    return tmp


def main():
    foldername = "grid_layout_map/"
    os.makedirs(foldername, exist_ok=True)
    foldername += str(time.time())
    # g = CombinMap()
    # g.addFace(8)
    # g.triangulate()
    colorsPerBiome = {
        "ile": 3,
        "plage": 4,
        "lagon": 2,
        "fonds": 1
    }
    sizeX, sizeY = 100, 100
    margin = 3
    rules_filename = "test_adjacency.json"
    model_filename = "test_model.json"
    with open(rules_filename, 'r') as rules_file, open(model_filename, 'r') as model_file:
        rules_content = json.load(rules_file)
        model_content = json.load(model_file)
        nodes = []
        nodesNames = []
        colors = []
        for model in model_content["subbiotopes"]:
            nodesToAdd = [model["model-name"] for _ in range(max(1, model["quantity"] + random.randint(1, 2)))]
            nodesNames += nodesToAdd
            nodes += nodesToAdd
            colors += [colorsPerBiome[n] for n in nodesToAdd]

        instancedNodes = []
        edge2node = []
        g = CombinMap()
        newEdges = set(g.addFace())
        for edge in newEdges:
            edge.affectSource(0)
        edge2node.append((newEdges, 0))
        # Add it 2 times, just there, so we can find 2 nodes' ID on each edges
        edge2node.append((newEdges, 0))
        instancedNodes.append(0)

        tries = 0
        totalIter = 0
        while len(instancedNodes) < len(nodes):
            totalIter += 1
            exteriorFace = g.exteriorFace()
            usedEdge = random.choice(exteriorFace)
            nodeA, nodeB, instanceID = getCompatibleNodesForNewFace(usedEdge, edge2node, instancedNodes, nodes, rules_content)
            if instanceID >= 0:
                tries = 0
                newEdges = g.addFace(adjacentBrins=[usedEdge])
                for e in newEdges:
                    if e.getAffectedDest() == -1:
                        e.affectDest(instanceID)
                        break  # The sources/destinations will repercute, no need for more checks
                edge2node.append((newEdges, instanceID))
                instancedNodes.append(instanceID)
            else:
                tries += 1
                if tries > 100:
                    tries = 0
                    # If too many tries have been done, duplicate all exterior edges, there should be at least
                    # one edge that works for new nodes
                    for e in exteriorFace:
                        previous, new = e.subdivide()
                        previous.affectDest(previous.getAffectedSource())
                    g.triangulate()
        doCollapse = True
        while doCollapse:
            doCollapse = False
            for brin in g.root.getAllBrins():
                # print(f"Do I collapse {brin}?", end = " ")
                if brin.getAffectedSource() == brin.getAffectedDest() and brin.beta1 != brin.beta2:
                    # print("Yes")
                    g.collapse(brin, mergeIfNeeded = False)
                    doCollapse = True
                    break

        # g.debug()
        grid = np.ones((sizeX, sizeY))
        nodes = g.affectNodes()
        edges = g.getCorrectionUnorientedEdges()
        targetSizes = [random.randint((sizeX * sizeY)//(len(nodes) + 5), (sizeX * sizeY)//(len(nodes))) for _ in range(len(nodes))]
        # targetSizes = [(sizeX * sizeY)//(len(nodes) + 2) for _ in range(len(nodes))]

        G, _ = g.toNetworkX()
        G_pos = nx.kamada_kawai_layout(G)

        # geom = g.toGeometricGraph()
        # geom.positions = geometryFromInstancingNodes(geom)
        # geom.positions = repositionNodesWithForces(geom)
        # G_pos = {i: geom.positions[i, :2] for i in range(len(geom.nodes))}

        minX, minY, maxX, maxY = math.inf, math.inf, -math.inf, -math.inf
        nodesPos = [[0, 0]] * len(nodes)
        for i, val in G_pos.items():
            nodesPos[i] = val.tolist()
            minX = min(minX, val[0])
            maxX = max(maxX, val[0])
            minY = min(minY, val[1])
            maxY = max(maxY, val[1])

        for node in nodesPos:
            node[0] = margin + int((sizeX-(margin * 2)) * (node[0] - minX) / (maxX - minX))
            node[1] = margin + int((sizeY-(margin * 2)) * (node[1] - minY) / (maxY - minY))

        # pfg = pathfinderGrid(matrix=grid)
        finder = BiAStarFinder(diagonal_movement=pathfinderDiagonalMovement.always)
        regionsSizes = getRegionsSizes(gridAsRegions(grid), len(nodes))
        previousRegionsSizes = [-1] * len(regionsSizes)
        tries = 100
        while np.any(grid[:] == 1) and \
                np.any(np.array(getRegionsSizes(gridAsRegions(grid), len(targetSizes))) < np.array(targetSizes)) and \
                np.any(np.array(regionsSizes) != np.array(previousRegionsSizes)) and tries > 0:
            tries -= 1
            keptEdges = []
            for i, (nodeA, nodeB) in enumerate(edges):
                tmpGrid = grid.copy()
                for x in range(sizeX):
                    for y in range(sizeY):
                        cell = tmpGrid[x, y]
                        if cell > 1:
                            if cell == nodeA or cell == nodeB:
                                tmpGrid[x, y] = 1
                            else:
                                tmpGrid[x, y] = -1

                pfg = pathfinderGrid(matrix = tmpGrid.T)
                start = pfg.node(*nodesPos[nodeA])
                end = pfg.node(*nodesPos[nodeB])
                path, runs = finder.find_path(start, end, pfg)
                if len(path) == 0:
                    keptEdges.append((nodeA, nodeB))
                pfg.cleanup()
                for cell in path[1:len(path)//2]:
                    grid[cell[0], cell[1]] = -(nodeA + 2)
                for cell in path[len(path)//2:-1]:
                    grid[cell[0], cell[1]] = -(nodeB + 2)
            grid = dilute(grid, targetSizes)
            previousRegionsSizes = regionsSizes
            regionsSizes = getRegionsSizes(gridAsRegions(grid), len(nodes))
            edges = keptEdges
        g.debug(labels= lambda i: str(nodesNames[i]) + "#" + str(i), block = False, tofile=foldername + "graph.png")
        # fig = plt.figure()
        gridShow(grid, colors, tofile=foldername + "result.png")


if __name__ == "__main__":
    # random.seed(420)
    for _ in range(100):
        main()
