import matplotlib.pyplot as plt
import matplotlib.colors
import os
import numpy as np
import pandas as pd
import karstnet as kn

nodeindex = 0


class Node:
    parent = None
    parents = None
    children = None
    pos = None
    index = 0
    displayable = False

    def __init__(self, _parent, _pos, _display):
        global nodeindex
        self.index = nodeindex + 1
        nodeindex += 1
        self.parent = _parent
        self.pos = []
        self.parents = []
        self.children = []
        self.displayable = _display
        for p in _pos:
            if not isinstance(p, float):
                p = float(p)
            self.pos.append(p)

def display_points(filename):
    nodes_filename = filename + "_nodes.dat"
    links_filename = filename + "_links.dat"
    nodes = []

    with open(nodes_filename, "r") as f:
        for line in f:
            splitted = line.split(" ")
            if len(splitted) > 0:
                nodes.append(Node(None, splitted[0:3], splitted[3] == "1"))

    with open(links_filename, "r") as f:
        for line in f:
            splitted = line.split(" ")
            # if int(splitted[1]) - 1 >= 18:
            #     nodes[int(splitted[1]) - 1].parent = nodes[int(splitted[0]) - 1]
            # else:
            nodes[int(splitted[1]) - 1].children.append(nodes[int(splitted[0]) - 1])
            nodes[int(splitted[0]) - 1].parents.append(nodes[int(splitted[1]) - 1])

    maxZ = -999999
    minZ = 999999
    for node in nodes:
        if maxZ < node.pos[2]:
            maxZ = node.pos[2]
        if minZ > node.pos[2]:
            minZ = node.pos[2]
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    i = 0
    nb_lines = 0
    for node in nodes:
        i += 1
        # color = matplotlib.colors.hsv_to_rgb((i / len(nodes), 1.0, 1.0))
        color = [1 - (node.pos[2] - minZ)/(maxZ - minZ), 0.0, (node.pos[2] - minZ)/(maxZ - minZ)]
        for parent in node.parents:
            nb_lines += 1
            ax.plot3D([node.pos[0], parent.pos[0]],
                      [node.pos[1], parent.pos[1]],
                      [node.pos[2], parent.pos[2]], c=color)
        if (node.displayable):
            print(node.index)
            ax.scatter3D([node.pos[0]], [node.pos[1]], [node.pos[2]], color=color, alpha=0.2)
            ax.text3D(node.pos[0], node.pos[1], node.pos[2], str(node.index-1), color=color)

    # print(nb_lines)
    plt.show()


def main():
    filename = "data/karst"
    """    # filename = "data/Huttes"
    # filename = "data/Sakany"
    # filename = "data/superimposed"
    # filename = "data/spongework"
    # filename = "data/rectilinear"
    # filename = "data/gorge"
    nodes_filename = filename + "_nodes.dat"
    links_filename = filename + "_links.dat"
    nodes = []

    nodes_dict = {}
    edges = []

    with open(nodes_filename, "r") as f:
        for line in f:
            splitted = line.split(" ")
            if len(splitted) == 3:
                nodes.append(Node(None, splitted, True))
                nodes_dict[nodes[-1].index] = nodes[-1].pos

    with open(links_filename, "r") as f:
        for line in f:
            splitted = line.split(" ")
            # if int(splitted[1]) - 1 >= 18:
            #     nodes[int(splitted[1]) - 1].parent = nodes[int(splitted[0]) - 1]
            # else:
            nodes[int(splitted[0]) - 1].parent = nodes[int(splitted[1]) - 1]
            edges.append((nodes[int(splitted[0]) - 1].index, nodes[int(splitted[0]) - 1].parent.index))
"""
    karst = kn.from_nodlink_dat(filename)
    results = karst.characterize_graph(verbose=True)
    karst.basic_analysis()
    karst.plot3()
    return
    """
    minX, maxX, minY, maxY, minZ, maxZ = None, None, None, None, None, None
    for node in nodes:
        minX = node.pos[0] if minX is None or minX > node.pos[0] else minX
        maxX = node.pos[0] if maxX is None or maxX < node.pos[0] else maxX
        minY = node.pos[1] if minY is None or minY > node.pos[1] else minY
        maxY = node.pos[1] if maxY is None or maxY < node.pos[1] else maxY
        minZ = node.pos[2] if minZ is None or minZ > node.pos[2] else minZ
        maxZ = node.pos[2] if maxZ is None or maxZ < node.pos[2] else maxZ

    print(minX, " - ", maxX)
    print(minY, " - ", maxY)
    print(minZ, " - ", maxZ)
    maxX -= minX
    maxY -= minY
    maxZ -= minZ
    maxY *= 8
    maxMax = 1.0  # max(maxX, max(maxY, maxZ))

    print("Box aspect : ", (maxX / maxMax, maxY / maxMax, maxZ / maxMax))
    # return
    fig = plt.figure()
    ax = plt.axes(projection="3d") #, box_aspect=(maxX / maxMax, maxY / maxMax, maxZ / maxMax))
    i = 0
    for node in nodes:
        i += 1
        color = matplotlib.colors.hsv_to_rgb((i / len(nodes), 1.0, 1.0))
        if node.parent is not None:
            ax.plot3D([node.pos[0], node.parent.pos[0]],
                      [node.pos[1], node.parent.pos[1]],
                      [node.pos[2], node.parent.pos[2]], c=color)
        else:
            print(node.index)
        # ax.scatter3D([node.pos[0]], [node.pos[1]], [node.pos[2]], color=color)
        # ax.text3D(node.pos[0], node.pos[1], node.pos[2], str(node.index), color=color)
    plt.show()"""


if __name__ == '__main__':
    main()
