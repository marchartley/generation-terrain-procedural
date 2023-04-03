from random import random
import math

radius = 10
width = 100
height = 100
depth = 100

def dist(vecA, vecB):
    return dist_2(vecA, vecB)**(1/2)

def dist_2(vecA, vecB):
    return (vecA[0] - vecB[0])**2 + (vecA[1] - vecB[1])**2 + (vecA[2] - vecB[2])**2

def make_poisson():
    nb_tries = 100
    current_tries = 0
    poisson_grid = []
    while True:
        if current_tries <= 0:
            break
        current_tries += 1
        p = [random() * width, random() * height, random() * depth]
        keep = True
        for other in poisson_grid:
            if dist_2(p, other) < radius**2:
                keep = False
                break
        if keep:
            poisson_grid.append(p)
            current_tries = nb_tries
    return poisson_grid

def make_voroinoi_as_grid():
    voronoi_grid = []
    for x in range(int(width/radius)):
        for y in range(int(height/radius)):
            for z in range(int(depth/radius)):
                voronoi_grid.append([
                    (random() + x) * radius,
                    (random() + y) * radius,
                    (random() + z) * radius
                ])
    return voronoi_grid

def make_voronoi():
    voronoi = []
    for i in range((width*height*depth)/(radius**3)):
        voronoi.append([random() * width, random() * height, random() * depth])

def compare_pathfinders():
    for nb_nodes in range(1, 100):
        nb_edges = 48*nb_nodes
        floyd = nb_nodes**3
        djikstra = nb_edges + nb_nodes * math.log2(nb_nodes)
        print("Pour", nb_nodes, "nodes,", ("Djikstra" if djikstra < floyd else "Floyd"), "est meilleur (D=", djikstra, "; F=", floyd, ")")


def main():
    compare_pathfinders()
    return
    poisson_distrib = make_poisson()
    voronoi_grid_distrib = make_voroinoi_as_grid()
    voronoi_distrib = make_voronoi()
    pass

if __name__ == "__main__":
    main()