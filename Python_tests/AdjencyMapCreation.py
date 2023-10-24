import math
import random
from typing import List, Tuple, Set

from GeometricGraph import GeometricGraph, geometryFromInstancingNodes, repositionNodesWithForces, \
    iterativeRepositionNodesWithForce
from RegionGrowing import Region
from TopoMap import CombinMap, Brin
import json
import networkx as nx
import matplotlib.pyplot as plt

from heightmap_extractor.constraints_tests import Vector2D

import os


def getNodesOnEdge(edge: Brin, edges2node: List[Tuple[Set[Brin], int]]) -> Tuple[int, int]:
    return edge.getAffectedSource(), edge.getAffectedDest()

def getCompatibleNodesForNewFace(edge: Brin,
                                 edges2node: List[Tuple[Set[Brin], int]],
                                 instancedNodes: List[int],
                                 allNodes: List[str],
                                 rules) -> Tuple[int, int, int]:
    nodeA, nodeB = getNodesOnEdge(edge, edges2node)
    # print(f"Instancing between {allNodes[nodeA]} and {allNodes[nodeB]}...", end=" ")
    if nodeA < 0 or nodeB < 0:
        # print("Error.")
        return nodeA, nodeB, -1
    availableNodes: List[int] = []
    withoutRules = []
    for i in range(len(allNodes)):
        if i not in instancedNodes:
            rule = rules[allNodes[i]]
            withoutRules.append(i)
            if allNodes[nodeA] in rule and allNodes[nodeB] in rule:
                availableNodes.append(i)
    # print(f"Without rule, {[allNodes[i] for i in withoutRules]} can be instanciated, "
    #       f"but {[allNodes[i] for i in availableNodes]} with them.")
    if len(availableNodes) == 0:
        return nodeA, nodeB, -1
    return nodeA, nodeB, random.choice(availableNodes)


def main():
    rules_filename = "test_adjacency.json"
    model_filename = "test_model.json"
    with open(rules_filename, 'r') as rules_file, open(model_filename, 'r') as model_file:
        rules_content = json.load(rules_file)
        model_content = json.load(model_file)

        first = True
        for _ in range(100):
            nodes = []
            for model in model_content["subbiotopes"]:
                nodes += [model["model-name"] for _ in range(max(1, model["quantity"] + random.randint(0, 8)))]

            # if first:
            #     first = False
            #     continue
            foldername = "./saved_runs/" + "".join([str(n[0].upper()) for n in nodes]) + "/"
            os.makedirs(foldername, exist_ok=True)

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
                    # else:
                    #     print("No")

            g.debug(lambda n: (str(nodes[n])[0].upper() + "#" + str(n)) if n < len(nodes) else "N", tofile=foldername + "PDF_0.pdf")
            g.debug(lambda n: (str(nodes[n])[0].upper() + "#" + str(n)) if n < len(nodes) else "N", tofile=foldername + "IMG_0.png")
            g.addNeutralComponents(nodes)
            g.debug(lambda n: (str(nodes[n])[0].upper() + "#" + str(n)) if n < len(nodes) else "N", tofile=foldername + "PDF_1.pdf")
            g.debug(lambda n: (str(nodes[n])[0].upper() + "#" + str(n)) if n < len(nodes) else "N", tofile=foldername + "IMG_1.png")
            for mainIter in range(3):
                graph = g.toGeometricGraph()
                graph.positions = geometryFromInstancingNodes(graph)
                previousPositions = graph.positions.copy()
                nbIterations = 100
                for iteration in range(nbIterations):
                    # plt.pause(0.05)
                    plt.clf()
                    graph.positions = iterativeRepositionNodesWithForce(graph, deltaMove=0.8)
                    diff = graph.positions - previousPositions
                    if np.sum(abs(diff)) / len(graph.nodes) < 1e-5:
                        break
                    previousPositions = graph.positions.copy()
                    graph.draw(withLabels=True, withVoronoi=True, withDelaunay=True, title=f"Step {iteration + 1}/{nbIterations}")
                    # plt.draw()
                    plt.savefig(foldername + f"PDF_main{mainIter}-iterA{iteration + 1}.pdf")
                    plt.savefig(foldername + f"IMG_main{mainIter}-iterA{iteration + 1}.png")

                networkXpositions = nx.kamada_kawai_layout(g.toNetworkX()[0])
                # Region growing step:
                randomAreas = [(sum([weight for n, weight in graph.nodes[i].connectedNodes]) / len(
                    graph.nodes[i].connectedNodes)) ** 2 * math.pi for i in range(len(graph.nodes))]
                # regions = [Region(Vector2D(graph.positions[i, 0], graph.positions[i, 1]), randomAreas[i], 30, i) for i in
                #            range(len(graph.nodes))]
                print(networkXpositions)
                regions = [Region(Vector2D(networkXpositions[i][0], networkXpositions[i][1]), randomAreas[i], 30, i) for i in
                           range(len(graph.nodes))]

                for i, region in enumerate(regions):
                    # region.initVerticesFromVoronoi(graph, i)
                    region.initVerticesFromNodeAndNeighbors(graph, i)
                    region.initVerticesAsUnitCircle(graph, i, 0.05)

                # for i, region in enumerate(regions):
                #     region.draw(withAreas=randomAreas[i], withVertices=False, withStrengths=False)
                # graph.draw(title="Initial graph", withLabels=False, withVoronoi=False)
                # plt.draw()

                for iteration in range(nbIterations):
                    for i, region in enumerate(regions):
                        if i < len(nodes):
                            region.setVerticesNormalsTowardOtherRegion(graph, i, regions)
                            region.setVerticesStrengthTowardNodes(graph, i, regions, 0.05)
                        else:
                            region.setVerticesStrength([1.0] * len(region.vertices), 0.03)

                    # plt.pause(0.01)
                    plt.clf()
                    allRegionFull = True
                    for i, region in enumerate(regions):
                        region.growWithNeighborRegionDisplacement(regions)
                        # region.smooth(0.5)
                        if i >= len(nodes):
                            continue
                        region.draw(withAreas=randomAreas[i], withVertices=False, withStrengths=False)
                        if region.getProportionAreaCovered() < 1.0:
                            allRegionFull = False
                    graph.draw(title=f"Initial graph (iter {iteration + 1}/{nbIterations})", withLabels=False,
                               withVoronoi=False, selectedNodes=list(range(len(nodes))))
                    # plt.draw()
                    plt.savefig(foldername + f"PDF_main{mainIter}-iterB{iteration + 1}.pdf")
                    plt.savefig(foldername + f"IMG_main{mainIter}-iterB{iteration + 1}.png")
                    if allRegionFull:
                        break
                # plt.ioff()
                # plt.show()
                plt.clf()


if __name__ == "__main__":
    import numpy as np

    seed = 12480
    random.seed(seed)
    np.random.seed(seed)
    main()
