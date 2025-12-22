import matplotlib.pyplot as plt


grid = [0 for _ in range(9*9)]

grid[:9] = [i for i in range(9)]

def indx(x, y): return x * 9 + y

def print_grid(g):
    for col in range(9):
        for row in range(9):
            print(f"{g[indx(row, col)]} ", end="")
            if  (row + 1) % 3 == 0 :
                print("|", end="")
        print("")
        if (col + 1) % 3 == 0:
            print("-" * 9)


print_grid(grid)