import math
from typing import List, Tuple, Union, Optional

from Graph import Node
from heightmap_extractor.constraints_tests import Vector2D
import random
from GeometricGraph import GeometricGraph, geometryFromInstancingNodes, repositionNodesWithForces, line_intersection, \
    compute_dist
import matplotlib.pyplot as plt
import noise
import matplotlib.patches
import numpy as np
from pathfinding.core.grid import Grid as pathfinderGrid
from pathfinding.core.diagonal_movement import DiagonalMovement as pathfinderDiagonalMovement
from pathfinding.finder.bi_a_star import BiAStarFinder as BiAStarFinder


class Region:
    vertices: List[Vector2D] = []
    verticesStrength: List[float]
    fixedVertices: List[bool]
    id: int
    targetArea: float
    center: Vector2D
    precomputed_normals: List[Optional[Vector2D]] = []

    def __init__(self, center: Vector2D, targetArea: float = None, nbVertices: int = 20, _id: int = None):
        octaves = 6
        persistence = 0.5
        lacunarity = 2.0
        scale = 1.0
        meanTarget = .2
        self.center = center
        self.targetArea = targetArea
        self.id = _id
        self.vertices = [center.copy()
                         + Vector2D(math.cos(i * 2 * math.pi / nbVertices),
                                    math.sin(i * 2 * math.pi / nbVertices)) * .01
                         for i in range(nbVertices)]
        self.verticesStrength = [noise.pnoise2(self.vertices[i].x * scale, self.vertices[i].y * scale,
                                               octaves=octaves,
                                               persistence=persistence,
                                               lacunarity=lacunarity,
                                               base=0) for i in range(nbVertices)]
        # _min, _max = min(self.verticesStrength), max(self.verticesStrength)
        # self.verticesStrength = [(vert - _min)/(_max - _min) for vert in self.verticesStrength]
        self.verticesStrength = [abs(s) for s in self.verticesStrength]
        self.verticesStrength = [1.0 for vert in self.verticesStrength]

        _mean = sum(self.verticesStrength) / len(self.verticesStrength)
        self.verticesStrength = [vert * (meanTarget / _mean) for vert in self.verticesStrength]
        self.fixedVertices = [False for _ in self.vertices]
        self.precomputed_normals = [None for _ in self.vertices]

    def setVerticesStrength(self, newStrengths: List[float], meanTarget: float):
        # _min, _max = min(self.verticesStrength), max(self.verticesStrength)
        # self.verticesStrength = [(vert - _min)/(_max - _min) for vert in self.verticesStrength]
        # self.verticesStrength = [abs(s) for s in self.verticesStrength]
        # self.verticesStrength = [1.0 for vert in self.verticesStrength]

        _mean = sum(newStrengths) / len(newStrengths)
        self.verticesStrength = [vert * (meanTarget / _mean) for vert in newStrengths]

    def setVerticesStrengthTowardNodes(self, graph: GeometricGraph, selfIndex: int, otherRegions: List['Region'], meanStrength: float = 1.0):
        newStrengths = [1.0 for _ in range(len(self.vertices))]
        consideredNeighbors = []
        for neighbor, weight in graph.nodes[selfIndex].connectedNodes:
            if not self.alreadyAdjacentTo(otherRegions[neighbor.id]):
                consideredNeighbors.append((neighbor.id, weight))
        if not consideredNeighbors:
            self.setVerticesStrength(newStrengths, meanStrength)
            return

        closestVertexPerRegion: List[List[Optional[Vector2D]]] = [[None for _ in range(len(self.vertices))] for _ in range(len(graph.nodes[selfIndex].connectedNodes))]
        furthestVertexPerRegion: List[List[Optional[Vector2D]]] = [[None for _ in range(len(self.vertices))] for _ in range(len(graph.nodes[selfIndex].connectedNodes))]
        for i, (neighbor, weight) in enumerate(graph.nodes[selfIndex].connectedNodes):
            if self.alreadyAdjacentTo(otherRegions[neighbor.id]):
                continue
            for iDir, vertex in enumerate(self.vertices):
                for oVertex in otherRegions[neighbor.id].vertices:
                    nodeDist = (oVertex - vertex)
                    if closestVertexPerRegion[i][iDir] is None or nodeDist.norm2() < (closestVertexPerRegion[i][iDir] - vertex).norm2():
                        closestVertexPerRegion[i][iDir] = oVertex
                    if furthestVertexPerRegion[i][iDir] is None or nodeDist.norm2() > (furthestVertexPerRegion[i][iDir] - vertex).norm2():
                        furthestVertexPerRegion[i][iDir] = oVertex

            for iDir, (closeVertex, farVertex) in enumerate(zip(closestVertexPerRegion[i], furthestVertexPerRegion[i])):
                vertex = self.vertices[iDir]
                if closeVertex is None or farVertex is None:
                    continue
                nodeDist = (farVertex - vertex) # ((closeVertex - vertex) + (farVertex - vertex)) / 2
                nodeDir = nodeDist.normalized()
                addStrength = (nodeDir.dot(self.normal(iDir))) #  * nodeDist.norm())
                newStrengths[iDir] += max(0.0, addStrength)
        self.setVerticesStrength(newStrengths, meanStrength)
        """        consideredNeighbors = []
        for neighbor, weight in graph.nodes[selfIndex].connectedNodes:
            if not self.alreadyAdjacentTo(otherRegions[neighbor.id]):
                consideredNeighbors.append((neighbor.id, weight))

        newStrengths = []
        for iDir, vertex in enumerate(self.vertices):
            direction = self.normal(iDir).normalized()
            strength = 1.0
            for neighbor, weight in consideredNeighbors:
                nodeDir = (self.center - Vector2D(graph.positions[neighbor, :2])).normalized()
                addStrength = (direction.dot(nodeDir) * weight)
                strength += max(0.0, addStrength)
            newStrengths.append(strength)
        self.setVerticesStrength(newStrengths, 0.5)"""

    def setVerticesNormalsTowardNodes(self, graph: GeometricGraph, selfIndex: int, otherRegions: List['Region']):
        newNormals = []
        for iDir, vertex in enumerate(self.vertices):
            direction = self.normal(iDir).normalized()
            newDir = direction.copy()
            for neighbor, weight in graph.nodes[selfIndex].connectedNodes:
                nodeDist = (otherRegions[neighbor.id].center - self.center)
                nodeDir = nodeDist.normalized()
                addStrength = (nodeDir.dot(direction) * nodeDist.norm2())
                newDir += direction * max(0.0, addStrength)
            newNormals.append(newDir.normalized())
        self.precomputed_normals = newNormals

    def setVerticesNormalsTowardOtherRegion(self, graph: GeometricGraph, selfIndex: int, otherRegions: List['Region'], distortionFactor: float = .5):
        newNormals = [Vector2D() for _ in range(len(self.vertices))]
        closestVertexPerRegion: List[List[Optional[Vector2D]]] = [[None for _ in range(len(self.vertices))] for _ in range(len(graph.nodes[selfIndex].connectedNodes))]
        furthestVertexPerRegion: List[List[Optional[Vector2D]]] = [[None for _ in range(len(self.vertices))] for _ in range(len(graph.nodes[selfIndex].connectedNodes))]
        for i, (neighbor, weight) in enumerate(graph.nodes[selfIndex].connectedNodes):
            if self.alreadyAdjacentTo(otherRegions[neighbor.id]):
                continue
            for iDir, vertex in enumerate(self.vertices):
                for oVertex in otherRegions[neighbor.id].vertices:
                    nodeDist = (oVertex - vertex)
                    if closestVertexPerRegion[i][iDir] is None or nodeDist.norm2() < (closestVertexPerRegion[i][iDir] - vertex).norm2():
                        closestVertexPerRegion[i][iDir] = oVertex
                    if furthestVertexPerRegion[i][iDir] is None or nodeDist.norm2() > (furthestVertexPerRegion[i][iDir] - vertex).norm2():
                        furthestVertexPerRegion[i][iDir] = oVertex

            for iDir, (closeVertex, farVertex) in enumerate(zip(closestVertexPerRegion[i], furthestVertexPerRegion[i])):
                vertex = self.vertices[iDir]
                if closeVertex is None or farVertex is None:
                    continue
                nodeDist = (farVertex - vertex) # ((closeVertex - vertex) + (farVertex - vertex)) / 2
                nodeDir = nodeDist.normalized()
                addStrength = (nodeDir.dot(self.normal(iDir))) #  * nodeDist.norm())
                newNormals[iDir] += nodeDir * max(0.0, addStrength)
        self.precomputed_normals = [(self.normal(i, False) * (1 - distortionFactor) + (n * distortionFactor).normalize()).normalize() for i, n in enumerate(newNormals)]

    def setVerticesNormalsTowardOtherRegionWithPathfinder(self, graph: GeometricGraph, selfIndex: int, otherRegions: List['Region'], distortionFactor: float = .5):
        newNormals = [Vector2D() for _ in range(len(self.vertices))]
        closestVertexPerRegion: List[List[Optional[Vector2D]]] = [[None for _ in range(len(self.vertices))] for _ in range(len(graph.nodes[selfIndex].connectedNodes))]
        furthestVertexPerRegion: List[List[Optional[Vector2D]]] = [[None for _ in range(len(self.vertices))] for _ in range(len(graph.nodes[selfIndex].connectedNodes))]
        for i, (neighbor, weight) in enumerate(graph.nodes[selfIndex].connectedNodes):
            if self.alreadyAdjacentTo(otherRegions[neighbor.id]):
                continue

            for iDir, vertex in enumerate(self.vertices):
                for oVertex in otherRegions[neighbor.id].vertices:
                    nodeDist = (oVertex - vertex)
                    if closestVertexPerRegion[i][iDir] is None or nodeDist.norm2() < (closestVertexPerRegion[i][iDir] - vertex).norm2():
                        closestVertexPerRegion[i][iDir] = oVertex
                    if furthestVertexPerRegion[i][iDir] is None or nodeDist.norm2() > (furthestVertexPerRegion[i][iDir] - vertex).norm2():
                        furthestVertexPerRegion[i][iDir] = oVertex

            dim = 20
            matrix = [[1 for _ in range(dim)] for _ in range(dim)]
            for x in range(dim):
                for y in range(dim):
                    for region in otherRegions:
                        marker = 0
                        if region.id is self.id:
                            marker = 1
                        elif region.id == neighbor.id:
                            marker = 2
                        if region.contains(self.center + Vector2D(x - dim/2, y - dim/2) * dim):
                            matrix[y][x] = marker

            grid = pathfinderGrid(matrix=matrix)
            finder = BiAStarFinder(diagonal_movement=pathfinderDiagonalMovement.only_when_no_obstacle)
            for iDir in range(len(newNormals)):
                startVec = (self.vertices[iDir] - self.center) + Vector2D(dim/2, dim/2)
                endVec = ((closestVertexPerRegion[i][iDir] - self.center) / 4) + Vector2D(dim/2, dim/2)
                s, e = grid.node(int(startVec.x), int(startVec.y)), grid.node(int(endVec.x), int(endVec.y))
                path, runs = finder.find_path(s, e, grid)
                grid.cleanup()
                if path:
                    midPath = path[len(path)//2]
                    newNormals[iDir] += Vector2D(midPath) - startVec
        self.precomputed_normals = [(self.normal(i, False) * (1 - distortionFactor) + (n * distortionFactor).normalize()).normalize() for i, n in enumerate(newNormals)]

    def alreadyAdjacentTo(self, otherRegion: 'Region', tol: float = 1.0):
        for iVert, vert in enumerate(self.vertices):
            proj = otherRegion.closestBorder(vert, evenIfAlreadyOutside=True)
            if (proj - vert).norm2() < tol*tol:
                return True
        return False

    def aabbox(self) -> Tuple[float, float, float, float]:
        x, y = [v.x for v in self.vertices], [v.y for v in self.vertices]
        # Return X, Y, W, H
        return min(x), min(y), max(x) - min(x), max(y) - min(y)

    def contains(self, point: Vector2D) -> bool:
        x, y, w, h = self.aabbox()
        if point.x < x or x + w < point.x or point.y < y or y + h < point.y:
            return False
        return pointInPolygon(point, self.vertices)

    def getIntersection(self, P0: Vector2D, P1: Vector2D) -> Tuple[int, Vector2D]:
        point = P0
        ray = P1 - P0
        intersectionPoints = []
        for i in range(len(self.vertices)):
            intersect = line_intersection(point, point + ray,
                                          self.vertices[i], self.vertices[(i + 1) % len(self.vertices)])
            if intersect:
                intersectionPoints.append((i, Vector2D(intersect[0], intersect[1])))

        closestDist = math.inf
        closestPoint = None
        closestIndex = None
        for i, intersect in intersectionPoints:
            if (P0 - intersect).norm2() < closestDist:
                closestDist = (P0 - intersect).norm2()
                closestPoint = intersect
                closestIndex = i
        return closestIndex, closestPoint

    def initVerticesAsUnitCircle(self, graph: GeometricGraph, nodeIndex: int, unit_size: float = 1.0):
        center = Vector2D(graph.positions[nodeIndex, :2])
        self.center = center
        for i, pos in enumerate(self.vertices):
            pos = Vector2D(math.cos(i * 2 * math.pi / len(self.vertices)),
                           math.sin(i * 2 * math.pi / len(self.vertices))) * unit_size
            self.vertices[i] = pos + center

    def initVerticesFromNodeAndNeighbors(self, graph: GeometricGraph, nodeIndex: int):
        center = Vector2D(graph.positions[nodeIndex, :2])
        self.center = center
        for i, pos in enumerate(self.vertices):
            pos = Vector2D(math.cos(i * 2 * math.pi / len(self.vertices)),
                           math.sin(i * 2 * math.pi / len(self.vertices)))
            # mult = 1.0
            # mostAlignedDist = 0.0
            mostAlignedNeighbor = None
            bestDot = 0
            for neighbor, weight in graph.nodes[nodeIndex].connectedNodes:
                neighborPos = Vector2D(graph.positions[neighbor.id, :2]) - center
                normalizedPos = neighborPos.normalized()
                dot = normalizedPos.dot(pos)
                if dot > bestDot:
                    bestDot = dot
                    mostAlignedNeighbor = neighborPos.copy()
            if bestDot > 0:
                pos *= bestDot * mostAlignedNeighbor.norm() * .1
            self.vertices[i] = pos + center

    def initVerticesFromVoronoi(self, graph: GeometricGraph, nodeIndex: int):
        voronoiResults = graph.computeVoronoi()
        if voronoiResults is not None:
            _, _, areas, voro = voronoiResults
            center = Vector2D(graph.positions[nodeIndex, :2])

            self.center = center
            for i, pos in enumerate(self.vertices):
                area = areas[nodeIndex]
                pos = Vector2D(math.cos(i * 2 * math.pi / len(self.vertices)),
                               math.sin(i * 2 * math.pi / len(self.vertices))) * 10 + center
                for iEdge in range(len(area)):
                    start, end = area[iEdge - 1], area[iEdge]
                    inter = line_intersection(center, pos,
                                              start, end)
                    if inter:
                        pos = Vector2D(inter)
                        # if len(voro.vertices) < max(posIdx):
                        #     pos = pos.normalized() * 2.0
                self.vertices[i] = pos

    def draw(self, withAreas: Union[bool, float] = False, withVertices: bool = True, withStrengths: bool = False):
        p = matplotlib.patches.Polygon([(v.x, v.y) for v in self.vertices], facecolor=(.5, .3, .3, .1))
        plt.gca().add_patch(p)
        line, = plt.plot([vec.x for vec in self.vertices + [self.vertices[0]]],
                         [vec.y for vec in self.vertices + [self.vertices[0]]])
        used_color = line.get_color()
        if withVertices:
            plt.scatter([vec.x for vec in self.vertices], [vec.y for vec in self.vertices], marker="x")
            for i in range(len(self.vertices)):
                plt.text(self.vertices[i].x, self.vertices[i].y, str(i))
        if withStrengths:
            maxSize = max(self.verticesStrength)
            targetSize = 1.50
            sizes = [targetSize * s / maxSize * (1 if self.fixedVertices[i] else 1) for i, s in enumerate(self.verticesStrength)]
            # plt.scatter([vec.x for vec in self.vertices], [vec.y for vec in self.vertices], sizes)
            normals = [self.normal(i) for i in range(len(self.vertices))]
            for i, normal in enumerate(normals):
                plt.plot([self.vertices[i].x, self.vertices[i].x + normals[i].x * sizes[i]],
                         [self.vertices[i].y, self.vertices[i].y + normals[i].y * sizes[i]],
                         "-", c=used_color)

        if withAreas:
            center = [sum([v.x for v in self.vertices]) / len(self.vertices),
                      sum([v.y for v in self.vertices]) / len(self.vertices)]
            if isinstance(withAreas, float):
                progress = self.area() / withAreas
                plt.text(center[0], center[1], str(round(self.area(), 2)) + "\n(" + str(int(progress * 100)) + "%)")
            else:
                plt.text(center[0], center[1], str(round(self.area(), 2)))

    def findSelfIntersections(self) -> List[Tuple[int, int, Vector2D]]:
        intersections = []
        for i in range(len(self.vertices)):
            i_plus_one = (i + 1) % len(self.vertices)
            for j in range(i + 1, len(self.vertices)):
                j_plus_one = (j + 1) % len(self.vertices)
                if i == j or i == j_plus_one or j == i_plus_one:
                    continue
                inter = line_intersection(self.vertices[i], self.vertices[i_plus_one],
                                          self.vertices[j], self.vertices[j_plus_one])
                if inter:
                    inter = Vector2D(inter[0], inter[1])
                    if j - i > len(self.vertices) / 2:
                        j -= len(self.vertices)
                        intersections.append((j + 1, i - 1, inter))
                    else:
                        intersections.append((i, j, inter))

        return intersections

    def findAndFilterSelfIntersections(self) -> List[Tuple[int, int, Vector2D]]:
        intersections = []
        allIntersections = self.findSelfIntersections()
        # allIntersections = list(filter(lambda tuple: tuple[1] - tuple[0] < len(self.vertices)/2, allIntersections))
        while allIntersections:
            firstIntersect = allIntersections[0]
            del allIntersections[0]
            start, end, vec = firstIntersect

            tmp = []
            for oStart, oEnd, oVec in allIntersections:
                if oStart > start and end > oEnd:
                    continue
                tmp.append((oStart, oEnd, oVec))
            allIntersections = tmp

            """# if (end - start) > len(self.vertices) / 2:
            #     continue #start = -start
            while allIntersections:
                otherIntersection = allIntersections[0]
                oStart, oEnd, _ = otherIntersection
                if oStart > start:
                    break
                else:
                    end = oEnd
                    del allIntersections[0]
            tmpAllIntersections = list(allIntersections)
            while tmpAllIntersections:
                otherIntersection = tmpAllIntersections[0]
                oStart, oEnd, _ = otherIntersection
                if oStart < end:
                    del tmpAllIntersections[0]
                elif oStart > end:
                    break
                else:
                    end = oEnd
                    del tmpAllIntersections[0]"""
            if end - start < len(self.vertices) / 3:
                intersections.append((start, end, vec))
        return intersections

    def cancelSelfIntersections(self):
        self.optimDone += 1
        if self.id == 3 and self.optimDone == 6:
            a = -1
        intersections = self.findAndFilterSelfIntersections()
        # intersections = self.findSelfIntersections()
        for iStart, iEnd, inter in intersections:
            for i in range(iStart, iEnd + 1):
                self.vertices[i % len(self.vertices)] = \
                    lerpOnPath((i - iStart) / ((iEnd + 1) - iStart),
                               [self.vertices[iStart - 1],
                                self.vertices[(iEnd + 2) % len(self.vertices)]])

    optimDone = 0

    def optimizeFrozenVertices(self):
        self.optimDone += 1
        freeRanges = []
        frozenRanges = []
        freeVert = []
        frozenVert = []
        for i in range(len(self.vertices)):
            if not self.fixedVertices[i]:
                freeVert.append(i)
                if len(freeRanges) > 0 and freeRanges[-1][1] == i - 1:
                    freeRanges[-1] = (freeRanges[-1][0], i)
                else:
                    freeRanges.append((i, i))
            else:
                frozenVert.append(i)
                if len(frozenRanges) > 0 and frozenRanges[-1][1] == i - 1:
                    if frozenRanges[-1][0] == frozenRanges[-1][1] \
                            or self.vertices[i].alignedWith(self.vertices[frozenRanges[-1][0]],
                                                            self.vertices[frozenRanges[-1][1]],
                                                            0.01):
                        frozenRanges[-1] = (frozenRanges[-1][0], i)
                    else:
                        frozenRanges.append((i, i))
                else:
                    frozenRanges.append((i, i))
        """
        if len(frozenRanges) > 1 and frozenRanges[-1][1] == len(self.vertices) - 1:
            newStart, newEnd = frozenRanges[-1][0], frozenRanges[0][1] - len(self.vertices)
            frozenRanges[-1] = (newStart, newEnd)
            del frozenRanges[0]
        if len(freeRanges) > 1 and freeRanges[-1][1] == len(self.vertices) - 1:
            newStart, newEnd = freeRanges[-1][0], freeRanges[0][1] - len(self.vertices)
            freeRanges[-1] = (newStart, newEnd)
            del freeRanges[0]"""

        frozenDuo = len([r for r in frozenRanges if r[0] != r[1]])
        frozenSolo = len([r for r in frozenRanges if r[0] == r[1]])
        frozenExtremes = (2 * frozenDuo + frozenSolo)

        freeDuo = len([r for r in freeRanges if r[0] != r[1]])
        freeSolo = len([r for r in freeRanges if r[0] == r[1]])
        freeExtremes = (2 * freeDuo + freeSolo)

        nbVerticesToReplace = len(self.vertices) - (frozenExtremes + freeSolo)
        nbVerticesPerFreeRange = ((nbVerticesToReplace / freeDuo) if freeDuo > 0 else 0)
        newVertices = []
        newFixed = []
        i = 0
        while len(freeRanges) or len(frozenRanges):
            if (freeRanges and frozenRanges and freeRanges[0][0] < frozenRanges[0][0]) or \
                    (freeRanges and not frozenRanges):
                start, end = freeRanges[0]
                if start == end:
                    newVertices.append(self.vertices[start].copy())
                    newFixed.append(False)
                    # nbVerticesToReplace -= 1
                else:
                    nbVerticesToInsert = int(nbVerticesPerFreeRange)
                    if 0 < nbVerticesToReplace - nbVerticesToInsert < nbVerticesToInsert:
                        nbVerticesToInsert = nbVerticesToReplace
                    path = [self.vertices[start - 1]] + [vec for vec in self.vertices[
                                                                        start:end + 1]]  # + [self.vertices[(end + 2) % len(self.vertices)]]
                    for i in range(nbVerticesToInsert):
                        newVertices.append(
                            lerpOnPath(i / ((nbVerticesToInsert - 1) if nbVerticesToInsert > 1 else 2.0), path))
                        newFixed.append(False)
                    nbVerticesToReplace -= nbVerticesToInsert
                del freeRanges[0]
            else:
                start, end = frozenRanges[0]
                if start == end:
                    newVertices.append(self.vertices[start])
                    newFixed.append(True)
                else:
                    newVertices.append(self.vertices[start])
                    newFixed.append(True)
                    newVertices.append(self.vertices[end])
                    newFixed.append(True)
                del frozenRanges[0]

        if self.id == 2:
            print(len(newVertices))
        self.vertices = newVertices
        self.fixedVertices = newFixed

    def closestBorder(self, point: Vector2D, normal: Vector2D = None, evenIfAlreadyOutside: bool = False) -> Vector2D:
        if not evenIfAlreadyOutside and not self.contains(point):
            return point
        closestPoint: Vector2D = Vector2D()
        closestDist = math.inf
        for i, vertex in enumerate(self.vertices):
            projectionPoint = None
            if normal is None:
                projectionPoint = projectPointOnLine(point, self.vertices[i - 1], self.vertices[i])
            else:
                inter = line_intersection(self.vertices[i - 1], self.vertices[i], point - normal * 100,
                                          point + normal * 100, 0.0001)
                if inter:
                    projectionPoint = Vector2D(inter[0], inter[1])
            if projectionPoint is None:
                continue
            dist = (point - projectionPoint).norm2()
            if dist < closestDist:
                closestPoint = projectionPoint.copy()
                closestDist = dist
        # if normal is not None and closestDist == math.inf:
        #     return self.closestBorder(point)
        return closestPoint

    def area(self) -> float:
        area = 0
        for i in range(len(self.vertices)):
            area += self.vertices[i].x * (
                    self.vertices[(i + 1) % len(self.vertices)].y - self.vertices[(i - 1) % len(self.vertices)].y)
        return abs(area) / 2.0

    def normal(self, vertIndex: int, usePrecomputedNormals: bool = True):
        if usePrecomputedNormals and len(self.precomputed_normals) > vertIndex and self.precomputed_normals[vertIndex] is not None:
            return self.precomputed_normals[vertIndex]
        previous = self.vertices[vertIndex - 1]
        nextVert = self.vertices[(vertIndex + 1) % len(self.vertices)]
        edge = (nextVert - previous).normalize()
        return edge.rotate(math.radians(-90))

    def getProportionAreaCovered(self) -> float:
        currentArea = self.area()
        return 0.0 if self.targetArea is None or self.targetArea == 0 else currentArea / self.targetArea

    def smooth(self, factor=1.0):
        for _ in range(math.ceil(factor)):
            newPositions = [vec.copy for vec in self.vertices]
            for i in range(len(newPositions)):
                if factor <= 1.0:
                    tmp = (self.vertices[i - 1] + self.vertices[(i + 1) % len(self.vertices)]) / 2.0
                    newPositions[i] = self.vertices[i] * (1.0 - factor) + tmp * factor
                else:
                    tmp = Vector2D()
                    for iVert in range(int(factor)):
                        tmp += (self.vertices[i - iVert] + self.vertices[(i + iVert) % len(self.vertices)])
                    newPositions[i] = tmp / (2 * int(factor))
            self.vertices = newPositions

    def growWithNeighborRegionDisplacement(self, otherRegions: List['Region']):
        self.fixedVertices = [False for _ in range(len(self.vertices))]
        if self.targetArea is not None and self.getProportionAreaCovered() >= 1.0:
            return
        newPositions = [vec.copy() for vec in self.vertices]
        for i, vert in enumerate(newPositions):
            if self.fixedVertices[i]:
                continue
            else:
                canMove = True
                direction = self.normal(i)
                predictedPos = vert + direction * self.verticesStrength[i]
                fixedPosition = None
                newVert = vert.copy()
                currentBest = math.inf
                for other in otherRegions:
                    if other is self:
                        continue
                    if other.getProportionAreaCovered() > self.getProportionAreaCovered():
                        continue
                    if other.contains(newVert):
                        test = other.closestBorder(newVert, direction * -1.0)
                        if (newPositions[i] - test).norm2() < currentBest:
                            currentBest = (newPositions[i] - test).norm2()
                            newVert = test.copy()
                if currentBest < math.inf:
                    newPositions[i] = newVert - direction * 0.001

                predictedPos = newPositions[i] + direction * self.verticesStrength[i]
                # Check if it will intersect another region
                for other in otherRegions:
                    if other is self:
                        continue
                    if other.contains(predictedPos):
                        idVertexIntersecting, possibleFixedPosition = other.getIntersection(vert, predictedPos)
                        ratio = (1.0 - (self.getProportionAreaCovered() / other.getProportionAreaCovered())
                                 if other.getProportionAreaCovered() > self.getProportionAreaCovered()
                                 else 0.0)
                        ratio = (1.0
                                 if other.getProportionAreaCovered() > self.getProportionAreaCovered()
                                 else 0.0)
                        if possibleFixedPosition is not None:
                            possibleFixedPosition = vert + (possibleFixedPosition - vert) * ratio

                        canMove = False
                        if possibleFixedPosition is not None:
                            if fixedPosition is None or (vert - possibleFixedPosition).norm2() < (
                                    vert - fixedPosition).norm2():
                                fixedPosition = possibleFixedPosition

                if self.contains(newPositions[i] + direction):
                    canMove = False

                if canMove:
                    newPositions[i] = predictedPos
                else:
                    if fixedPosition is not None:
                        newPositions[i] = fixedPosition
                    self.fixedVertices[i] = True

                for other in otherRegions:
                    if other is self or other.getProportionAreaCovered() < self.getProportionAreaCovered():
                        continue
                    for iVert, otherVert in enumerate(other.vertices):
                        if self.contains(otherVert):
                            other.vertices[iVert] = self.closestBorder(otherVert, other.normal(iVert) * -1.0)
        self.vertices = newPositions

    def growWithoutCheckingForOverlaps(self, otherRegions: List['Region']):
        self.fixedVertices = [False for _ in range(len(self.vertices))]
        if self.targetArea is not None and self.getProportionAreaCovered() >= 1.0:
            return
        newPositions = [vec.copy() for vec in self.vertices]
        for i, vert in enumerate(newPositions):
            if self.fixedVertices[i]:
                continue
            else:
                canMove = True
                direction = self.normal(i)
                fixedPosition = None
                newVert = vert.copy()
                currentBest = math.inf
                """for other in otherRegions:
                    if other is self:
                        continue
                    if other.getProportionAreaCovered() > self.getProportionAreaCovered():
                        continue
                    if other.contains(newVert):
                        test = other.closestBorder(newVert, direction * -1.0)
                        if (newPositions[i] - test).norm2() < currentBest:
                            currentBest = (newPositions[i] - test).norm2()
                            newVert = test.copy()
                if currentBest < math.inf:
                    newPositions[i] = newVert - direction * 0.001"""

                predictedPos = newPositions[i] + direction * self.verticesStrength[i]
                # Check if it will intersect another region
                """for other in otherRegions:
                    if other is self:
                        continue
                    if other.contains(predictedPos):
                        idVertexIntersecting, possibleFixedPosition = other.getIntersection(vert, predictedPos)
                        ratio = (1.0 - (self.getProportionAreaCovered() / other.getProportionAreaCovered())
                                 if other.getProportionAreaCovered() > self.getProportionAreaCovered()
                                 else 0.0)
                        ratio = (1.0
                                 if other.getProportionAreaCovered() > self.getProportionAreaCovered()
                                 else 0.0)
                        if possibleFixedPosition is not None:
                            possibleFixedPosition = vert + (possibleFixedPosition - vert) * ratio

                        # TODO : Check why I have better results without this line:
                        # canMove = False
                        if possibleFixedPosition is not None:
                            if fixedPosition is None or (vert - possibleFixedPosition).norm2() < (
                                    vert - fixedPosition).norm2():
                                fixedPosition = possibleFixedPosition"""

                if self.contains(newPositions[i] + direction):
                    canMove = False

                if canMove:
                    newPositions[i] = predictedPos
                else:
                    if fixedPosition is not None:
                        newPositions[i] = fixedPosition
                    self.fixedVertices[i] = True

                for other in otherRegions:
                    if other is self or other.getProportionAreaCovered() < self.getProportionAreaCovered():
                        continue
                    for iVert, otherVert in enumerate(other.vertices):
                        if self.contains(otherVert):
                            other.vertices[iVert] = self.closestBorder(otherVert, other.normal(iVert) * -1.0)
        self.vertices = newPositions

    def growNaive(self, otherRegions: List['Region']):
        if self.targetArea is not None and self.getProportionAreaCovered() >= 1.0:
            return
        newPositions = [vec.copy() for vec in self.vertices]
        if self.id == 1 and self.optimDone == 3:
            a = -1
        for i, vert in enumerate(newPositions):
            canMove = True
            direction = self.normal(i)
            predictedPos = vert + direction * self.verticesStrength[i]
            fixedPosition = None
            newVert = vert.copy()
            currentBest = math.inf
            for other in otherRegions:
                if other is self:
                    continue
                if other.contains(newVert):
                    test = other.closestBorder(newVert, direction * -1.0)
                    if (newPositions[i] - test).norm2() < currentBest:
                        currentBest = (newPositions[i] - test).norm2()
                        newVert = test.copy()
            if currentBest < math.inf:
                newPositions[i] = newVert - direction * 0.001

            predictedPos = newPositions[i] + direction * self.verticesStrength[i]
            # Check if it will intersect another region
            for other in otherRegions:
                if other is self:
                    continue
                if other.contains(predictedPos):
                    idVertexIntersecting, possibleFixedPosition = other.getIntersection(vert, predictedPos)
                    canMove = False
                    if possibleFixedPosition is not None:
                        if fixedPosition is None or (vert - possibleFixedPosition).norm2() < (
                                vert - fixedPosition).norm2():
                            fixedPosition = possibleFixedPosition.copy()

            if self.contains(newPositions[i] + direction):
                canMove = False

            if canMove:
                newPositions[i] = predictedPos
            else:
                if fixedPosition is not None:
                    newPositions[i] = fixedPosition
        self.vertices = newPositions


def projectPointOnLine(point: Vector2D, start: Vector2D, end: Vector2D, returnNoneIfOutside: bool = True) -> Optional[
    Vector2D]:
    AB = end - start
    t = AB.normalized().dot((point - start).normalize())
    if returnNoneIfOutside and (t < 0 or 1 < t):
        return None
    return start + AB * t


def pointInPolygon(point: Vector2D, polygon: List[Vector2D]) -> bool:
    ray = Vector2D(10000, 10)
    nbIntersections = 0
    for i in range(len(polygon)):
        if line_intersection(point, point + ray, polygon[i], polygon[(i + 1) % len(polygon)]):
            nbIntersections += 1
    return nbIntersections % 2 == 1


def lerpOnPath(x: float, path: List[Vector2D]) -> Vector2D:
    nbVert = len(path)
    X = x * (nbVert - 1)
    Xfloor = int(X)
    Xdec = X - Xfloor
    if x == 1.0:
        return path[-1].copy()
    return path[Xfloor] + (path[Xfloor + 1] - path[Xfloor]) * Xdec


def main():
    graph = GeometricGraph()
    while len(graph.nodes) == 0 or not graph.isPerfectlyPositioned():
        graph = GeometricGraph.randomPlanarGraph(8)
        graph.positions = geometryFromInstancingNodes(graph, 1.5, False)
        graph.positions = repositionNodesWithForces(graph)

    randomAreas = [(sum([weight for n, weight in graph.nodes[i].connectedNodes]) / len(
        graph.nodes[i].connectedNodes)) ** 2 * math.pi for i in range(len(graph.nodes))]
    regions = [Region(Vector2D(graph.positions[i, 0], graph.positions[i, 1]), randomAreas[i], 30, i) for i in
               range(len(graph.nodes))]

    for i, region in enumerate(regions):
        region.initVerticesFromVoronoi(graph, i)
        # region.initVerticesFromNodeAndNeighbors(graph, i)

    for i, region in enumerate(regions):
        region.draw(withAreas=randomAreas[i], withVertices=False, withStrengths=False)
    graph.draw(title="Initial graph", withLabels=False, withVoronoi=False)
    plt.draw()

    plt.ion()
    nbIterations = 100
    for iteration in range(nbIterations):
        for i, region in enumerate(regions):
            region.setVerticesNormalsTowardOtherRegion(graph, i, regions)
            region.setVerticesStrengthTowardNodes(graph, i, regions)

        plt.pause(0.01)
        plt.clf()
        for i, region in enumerate(regions):
            region.growWithoutCheckingForOverlaps(regions)
            region.smooth(0.5)
            region.draw(withAreas=randomAreas[i], withVertices=False, withStrengths=False)
        graph.draw(title=f"Initial graph (iter {iteration + 1}/{nbIterations})", withLabels=False, withVoronoi=False)
        plt.draw()
    plt.ioff()
    plt.show()


if __name__ == "__main__":
    random.seed(43)
    main()
