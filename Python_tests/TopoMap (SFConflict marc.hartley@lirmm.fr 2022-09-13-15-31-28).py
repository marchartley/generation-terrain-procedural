import math

from GeometricGraph import GeometricGraph
from Graph import Node, Graph, Vector2D
from typing import List, Union, Optional, Tuple, Set, Dict, Callable
import networkx as nx
import matplotlib.pyplot as plt

class Brin:
    currentBrinIndex: int = -1
    beta1: 'Brin' = None
    beta2: 'Brin' = None
    index: int = -1

    affectedSource: Optional[int] = -1
    affectedDest: Optional[int] = -1

    def __init__(self, index: int = -1):
        # self.beta1 = None
        # self.beta2 = None
        if index == -1:
            index = Brin.currentBrinIndex + 1
        if index > Brin.currentBrinIndex:
            Brin.currentBrinIndex = index
        self.index = index

    def getFace(self) -> List['Brin']:
        listOfBrins: List[Brin] = [self]
        current = self.beta1
        while current != self and current not in listOfBrins:
            listOfBrins.append(current)
            if current == current.beta1:
                return listOfBrins
            current = current.beta1
        return listOfBrins

    def previousBrin(self) -> 'Brin':
        return self.getFace()[-1]  # Return the last brin found in the face (knowing that the face starts with self)

    def getOppositeFace(self) -> List['Brin']:
        return self.beta2.getFace()

    def getTwoFaces(self) -> Tuple[List['Brin'], List['Brin']]:
        return self.getFace(), self.getOppositeFace()

    def getAllFaces(self) -> List[List['Brin']]:
        faces: List[List[Brin]] = []
        unseenBrins = self.getAllBrins()
        while len(unseenBrins) > 0:
            current = unseenBrins.pop()
            face = current.getFace()
            for brin in face:
                if brin in unseenBrins:
                    unseenBrins.remove(brin)
            if len(face) > 1:
                faces.append(face)
        return faces

    def getAllBrins(self) -> List['Brin']:
        seenBrins: List[Brin] = []
        unseenBrins: Set[Brin] = set()
        unseenBrins.add(self)

        while len(unseenBrins) > 0:
            current = unseenBrins.pop()
            seenBrins.append(current)
            if current.beta2 not in seenBrins:
                unseenBrins.add(current.beta2)
            if current.beta1 not in seenBrins:
                unseenBrins.add(current.beta1)

        return seenBrins

    def subdivide(self) -> Tuple['Brin', 'Brin']:
        nextBrin = self.beta1
        twinBrin = self.beta2
        # nextTwinBrin = twinBrin.beta1
        previousBrinTwin = twinBrin.previousBrin()

        self.beta1 = Brin()
        self.beta1.beta1 = nextBrin

        previousBrinTwin.beta1 = Brin()
        previousBrinTwin.beta1.beta1 = twinBrin

        previousBrinTwin.beta1.beta2 = self.beta1
        self.beta1.beta2 = previousBrinTwin.beta1

        return self, self.beta1

    def orbit(self) -> List['Brin']:
        if self.beta1 == self.beta1.beta1:
            return [self]
        brins: List[Brin] = [self]
        current = self.beta2.beta1
        while current != self:
            brins.append(current)
            if not current.isValid() or not current.beta2.isValid():
                return []
            else: # if current.beta1 != current.beta2:
                current = current.beta2.beta1
            # else:
            #     return []
        return brins

    def affectSource(self, node: int):
        for brin in self.orbit():
            # print(f"{brin} modified to {node}")
            brin.beta2.affectedDest = node
            brin.affectedSource = node

    def getAffectedSource(self) -> int:
        if self.affectedSource <= -1:
            # print(f"For {self}, checking {[self.beta2.orbit()]}")
            for brin in self.orbit():
                if brin.affectedSource > -1:
                    self.affectedSource = brin.affectedSource
                if brin.beta2.affectedDest > -1:
                    self.affectedSource = brin.beta2.affectedDest
        return self.affectedSource

    def affectDest(self, node: int):
        for brin in self.beta2.orbit():
            # print(f"{brin} modified to {node}")
            brin.affectedSource = node
            brin.beta2.affectedDest = node

    def getAffectedDest(self) -> int:
        if self.affectedDest <= -1:
            for brin in self.beta2.orbit():
                if brin.beta2.affectedDest > -1:
                    self.affectedDest = brin.beta2.affectedDest
                if brin.affectedSource > -1:
                    self.affectedDest = brin.affectedSource
        return self.affectedDest

    def affectSourceAndDest(self, source: int, dest: int):
        self.affectSource(source)
        self.affectDest(dest)

    def isValid(self) -> bool:
        return self.beta1 is not None and self.beta2 is not None

    def degree(self):
        return len(self.beta2.orbit())

    def __repr__(self):
        valid = self.isValid()
        return f"Brin #{self.index}, " \
               f"twin = {self.beta2.index if self.beta2 else 'none'}, " \
               f"next = {self.beta1.index if self.beta1 else 'none'} " \
               f"({self.getAffectedSource() if valid else '--'} -> {self.getAffectedDest() if valid else '--'})"

class CombinMap:
    root: Optional[Brin]

    def __init__(self):
        self.root = None

    def addFace(self, numberOfEdges: int = 3, adjacentBrins: List[Brin] = None) -> List[Brin]:
        if numberOfEdges < 3:
            raise Exception("We need at least 3 edges to create a new face.")
        if adjacentBrins is None:
            adjacentBrins = []
        nbEdgesToCreate = numberOfEdges - len(adjacentBrins)
        newBrins = [Brin() for _ in range(2 * nbEdgesToCreate)]

        if self.root is None:
            entryBrin = newBrins[0]
            exitBrin = newBrins[-1]
            self.root = newBrins[1]  # ID = 1 so it's outside, easier after, I guess
        elif len(adjacentBrins) > 0:
            entryBrin = adjacentBrins[0].previousBrin()
            exitBrin = adjacentBrins[-1].beta1
        else:
            raise Exception("Either you are creating a new graph, either you define at least one edge to paste the "
                            "new face to. ")
        # Twins and successions
        for i in range(len(newBrins) // 2):
            newBrins[2 * i].beta2 = newBrins[2 * i + 1]
            newBrins[2 * i + 1].beta2 = newBrins[2 * i]
            if i == 0:
                newBrins[1].beta1 = exitBrin
                # exitBrin.beta1 = newBrins[0]
                newBrins[0].beta1 = newBrins[2] if len(newBrins) > 2 else adjacentBrins[0]
            elif i == (len(newBrins) // 2 - 1):
                newBrins[2 * i].beta1 = entryBrin if len(adjacentBrins) == 0 else adjacentBrins[0]
                # entryBrin.beta1 = newBrins[2 * i + 1]
                newBrins[2 * i + 1].beta1 = newBrins[2 * i - 1]  # Firs condition imply that i >= 1
            else:
                newBrins[2 * i].beta1 = newBrins[2 * (i + 1)]
                newBrins[2 * i + 1].beta1 = newBrins[2 * i - 1]

        if len(adjacentBrins) > 0:
            adjacentBrins[-1].beta1 = newBrins[0]
            entryBrin.beta1 = newBrins[-1]
        return newBrins

    def allFaces(self) -> List[List[Brin]]:
        if self.root is not None:
            return self.root.getAllFaces()
        return []

    def exteriorFace(self) -> List[Brin]:
        faces = sorted(self.allFaces(), key = lambda l: len(l))
        return faces[-1]

    def triangulate(self):
        faces = sorted([f for f in self.allFaces() if len(f) > 3], key = lambda f: len(f))
        faces = faces[:-1]  # Do not triangulate the exterior face
        for face in faces:
            addedEdge = face[0]
            for i in range(2, len(face) - 1):
                newEdges = self.addFace(3, [addedEdge, face[i - 1]])
                addedEdge = newEdges[0] if newEdges[0] in addedEdge.beta2.orbit() else newEdges[1]

    def __repr__(self):
        return "Empty graph" if self.root is None else "Graph :\n" + "\n".join([str(b) for b in sorted(self.root.getAllBrins(), key = lambda brin: brin.index)])

    def getUnorientedEdges(self, ignoreIndices: int = math.inf) -> List[Tuple[int, int]]:
        allEdges = self.root.getAllBrins()
        allPreviousEdgesSource = [e.getAffectedSource() for e in allEdges]
        self.affectNodes()
        allNewEdgesSource = [e.getAffectedSource() for e in allEdges]
        unorientedEdges: List[Tuple[int, int]] = []
        unseenBrins = allEdges
        while len(unseenBrins) > 0:
            current = unseenBrins.pop()
            face = current.getFace()
            for brin in face:
                if brin in unseenBrins:
                    unseenBrins.remove(brin)
                if brin.beta2 in unseenBrins:
                    unseenBrins.remove(brin.beta2)
            for edge in face:
                unorientedEdges.append((edge.getAffectedSource(), edge.getAffectedDest()))

        # allEdges = self.root.getAllBrins()
        # for i, e in enumerate(allEdges):
        #     e.affectSource(allPreviousEdgesSource[i])
        return unorientedEdges

    def affectNodes(self) -> :
        allEdges = self.root.getAllBrins()
        nodeIndex = max([e.getAffectedSource() for e in allEdges] + [e.getAffectedDest() for e in allEdges])
        # nodeIndex = -1
        for e in allEdges:
            if e.getAffectedSource() == -1:
                nodeIndex += 1
                e.affectSource(nodeIndex)

    def toNetworkX(self) -> Tuple[nx.Graph, Dict[int, int]]:
        allEdges = self.root.getAllBrins()
        allPreviousEdgesSource = [e.getAffectedSource() for e in allEdges]
        allNewEdgesSource = [e.getAffectedSource() for e in allEdges]
        G = nx.Graph()
        G.add_edges_from(self.getUnorientedEdges())
        relabeling = {node: value for node, value in zip(allNewEdgesSource, allPreviousEdgesSource)}
        allEdges = self.root.getAllBrins()
        for i, e in enumerate(allEdges):
            e.affectSource(allPreviousEdgesSource[i])
        return G, relabeling

    def toGraph(self, ignoreIndices: int = math.inf) -> Graph:
        unorientedEdges = self.getUnorientedEdges(ignoreIndices)
        nodesID = list(range(1 + max([e[0] for e in unorientedEdges] + [e[1] for e in unorientedEdges]))) #list(set([e[0] for e in unorientedEdges] + [e[1] for e in unorientedEdges]))
        Node.currentIDs = []
        G = Graph([Node(_id) for _id in nodesID])
        for e in unorientedEdges:
            G.addNeighbors(e[0], e[1])
        for node in G.nodes:
            if len(node.connectedNodes) == 0:
                G.removeNode(node)
        return G

    def toGeometricGraph(self, ignoreIndices: int = math.inf) -> GeometricGraph:
        normalGraph = self.toGraph(ignoreIndices = ignoreIndices)
        return GeometricGraph.fromGraph(normalGraph)

    def debug(self, labels: Union[None, List, Callable, bool] = None, tofile: Union[bool, str] = False):
        G, _labels = self.toNetworkX()
        positions_kamada = nx.kamada_kawai_layout(G)
        usedLabels = {i: labels(l) for i, l in _labels.items()} if callable(labels) else labels
        nx.draw(G, positions_kamada, with_labels=False if isinstance(labels, bool) and not labels else True, labels = usedLabels)
        if not tofile:
            plt.show()
        else:
            plt.savefig(tofile)

    def clean(self):

        seenBrins: List[Brin] = []
        unseenBrins: Set[Brin] = set()
        unseenBrins.add(self.root)

        while len(unseenBrins) > 0:
            current = unseenBrins.pop()
            seenBrins.append(current)
            if current == current.beta1.beta1:
                twinA = current.beta2
                twinB = current.beta1.beta2
                if self.root in [current.beta1, current]:
                    self.root = twinA
                twinA.beta2 = twinB
                twinB.beta2 = twinA
            if current.beta2 not in seenBrins:
                unseenBrins.add(current.beta2)
            if current.beta1 not in seenBrins:
                unseenBrins.add(current.beta1)

    def collapse(self, brin: Brin, mergeIfNeeded: bool = False):
        brin.affectDest(brin.getAffectedSource())
        previousA = brin.previousBrin()
        previousB = brin.beta2.previousBrin()
        previousA.beta1 = brin.beta1
        previousB.beta1 = brin.beta2.beta1

        if self.root == brin or self.root == brin.beta2:
            self.root = previousA
        self.clean()
        """faces = brin.getAllFaces()
        exterior = sorted(faces, key=lambda face: len(face))[-1]
        if brin not in exterior and len(brin.getFace()) > 3:
            raise Exception()"""
        """twin = brin.beta2
        previousBrin = brin.previousBrin()
        previousTwin = twin.previousBrin()

        previousBrin.beta1 = brin.beta1
        previousTwin.beta1 = twin.beta1

        if brin == self.root or twin == self.root:
            self.root = self.root.beta1
        if mergeIfNeeded and previousBrin.beta1.beta1 is previousBrin:
            self.merge(previousBrin, previousBrin.beta1)"""

    # def createNeutralComponents(self):
    def merge(self, brinA: Brin, brinB: Brin):
        twinA = brinA.beta2
        twinB = brinB.beta2

        incomingA = twinA.previousBrin()
        incomingB = twinB.previousBrin()

        outgoingA = twinA.beta1
        outgoingB = twinB.beta1

        brinA.beta2 = brinB
        brinB.beta2 = brinA

        if incomingA != outgoingB:
            incomingA.beta1 = outgoingB
        else:
            self.collapse(incomingA) # incomingA.beta1 = outgoingB
        if incomingB != outgoingA:
            incomingB.beta1 = outgoingA
        else:
            self.collapse(incomingB)

        if self.root == twinA or self.root == twinB:
            self.root = brinA

    def addNeutralComponents(self, nodes) -> int:
        maxIndex = max([e.getAffectedSource() for e in self.root.getAllBrins()])
        neutralIndexLimit = maxIndex + 1
        while min([e.degree() for e in self.root.getAllBrins()]) <= 2:
            # self.debug(lambda n: (str(nodes[n])[0].upper() + "#" + str(n)) if n < len(nodes) else "N")

            # addNeutralBehind = False
            edges = self.exteriorFace()
            for i in range(len(edges)):
                while edges[i].degree() <= 2:
                    if edges[i].getAffectedSource() < neutralIndexLimit:
                        newEdges = self.addFace(3, [edges[i].previousBrin()])
                        for e in newEdges:
                            # Just catch the one which is the new "previousBrin()", change the source node and connect it
                            if e.beta1 == edges[i]:
                                maxIndex += 1
                                e.affectSource(maxIndex)
                                self.addFace(3, [e, edges[i]])
                                break
                    else:
                        # addNeutralBehind = True
                        newEdges = self.addFace(3, [edges[i].beta1])
                        for e in newEdges:
                            # Just catch the one which is the new "beta1", change the destination node and connect it
                            if edges[i].beta1 == e:
                                maxIndex += 1
                                e.affectDest(maxIndex)
                                self.addFace(3, [edges[i], e])
                                break
            # self.triangulate()

            self.clean()
            # self.debug(lambda n: (str(nodes[n])[0].upper() + "#" + str(n)) if n < len(nodes) else "N")

            # newExterior = self.exteriorFace()
            # for i in range(len(newExterior)):
            #     if newExterior[i].getAffectedSource() >= neutralIndexLimit and \
            #             newExterior[i].beta1.getAffectedDest() >= neutralIndexLimit:
            #         connect = self.addFace(3, [newExterior[i], newExterior[i].beta1])
            #         print(f"Collapsing {connect[0]}")
            #         self.collapse(connect[0])
            #         # self.merge(newExterior[i], newExterior[i].beta1)

            self.clean()
            continueCollapses = True
            while continueCollapses:
                self.clean()
                continueCollapses = False
                newExterior = self.exteriorFace()
                for i in range(len(newExterior)):
                    if newExterior[i].getAffectedSource() >= neutralIndexLimit and \
                            newExterior[i].getAffectedDest() >= neutralIndexLimit:
                        self.collapse(newExterior[i])
                        continueCollapses = True
                        break
            # self.debug(lambda n: (str(nodes[n])[0].upper() + "#" + str(n)) if n < len(nodes) else "N")

            self.clean()
            # if max([e.getAffectedSource() for e in self.root.getAllBrins()]) >= maxIndex:
            #     break
        return neutralIndexLimit


def main():
    g = CombinMap()
    g.addFace(4)
    g.affectNodes()
    # g.collapse(g.root)
    # g.collapse(g.root)
    g.debug(lambda l: l + 100)


    g = CombinMap()
    g.addFace(3)
    g.addFace(3, [g.exteriorFace()[0]])
    g.addFace(3, [g.exteriorFace()[0]])
    g.addFace(3, [g.exteriorFace()[0]])
    g.addFace(3, g.exteriorFace()[:2])
    for e in g.root.getAllBrins():
        e.affectSource(0)
    new = g.addFace(3, [g.exteriorFace()[0]])
    for e in new:
        if e.getAffectedSource() != 0:
            e.affectSource(1)
            break

    doCollapse = True
    while doCollapse:
        doCollapse = False
        for brin in g.root.getAllBrins():
            print(f"Do I collapse {brin}?")
            G, labels = g.toNetworkX()
            positions_kamada = nx.kamada_kawai_layout(G)
            positions_planar = nx.planar_layout(G)
            positions_spring = nx.spring_layout(G)
            positions_circular = nx.circular_layout(G)
            newLabels = {i: n for i, n in labels.items()}
            nx.draw(G, positions_kamada, labels=newLabels)  # , node_color=newColors))
            plt.show()
            if brin.getAffectedSource() == brin.getAffectedDest():
                g.collapse(brin, mergeIfNeeded=True)
                doCollapse = True
                break

    G, labels = g.toNetworkX()
    positions_kamada = nx.kamada_kawai_layout(G)
    positions_planar = nx.planar_layout(G)
    positions_spring = nx.spring_layout(G)
    positions_circular = nx.circular_layout(G)
    newLabels = {i: n for i, n in labels.items()}
    nx.draw(G, positions_kamada, labels=newLabels)  # , node_color=newColors))
    plt.show()
    # g.root.beta2.subdivide()
    # g.triangulate()
    # g.collapse(g.root)
    # g.collapse(g.root.beta2)
    # g.merge(g.root, g.root.beta1.beta1)
    faces = g.allFaces()
    # for i, face in enumerate(faces):
    #     print(f"Face #{i + 1} :\n" + "\n".join([str(b) for b in face]) + "\n")

    g.root.affectSource(42)
    # g.root.affectDest(0)
    # g.root.beta2.affectNode(0)
    # for edge in g.root.getAllBrins():
    #     print(edge, edge.getAffectedSource(), edge.getAffectedDest())
    G, _ = g.toNetworkX()
    nx.draw_planar(G, with_labels = True)
    plt.show()
    for edge in g.root.getAllBrins():
        print(edge, edge.getAffectedSource(), edge.getAffectedDest())
    print(g)


if __name__ == "__main__":
    main()