import matplotlib.pyplot as plt


def main():
    # name,     porosity,   grain-size,     density
    materials = [
        ["air",  1.00, 0.0001, 1.0],
        ["water",  1.00, 1.000, 1000.0],
        ["clay", 0.55, 0.0002, 1200.0],
        ["silt", 0.45, 0.0250, 1350.0],
        ["fine sand", 0.44, 0.2000, 1400.0],
        ["coarse sand", 0.43, 0.800, 1200],
        ["sand", 0.40, 1.000, 1200],
        ["gravel", 0.30, 2.000, 2000.0],
        ["rocks", 0.10, 0.000, 2100.0],
        ["limestone", 0.08, 0.0, 2700.0],
        ["sandstone", 0.20, 1.000, 2300.0],
        ["dolomite", 0.05, 0.20, 2700.0],
        ["shale", 0.15, 0.004, 2400.0],
    ]
    x = [m[1] for m in materials]
    y = [m[2] for m in materials]
    z = [m[3] for m in materials]
    s = [m[0] for m in materials]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z)
    ax.set_xlabel("Porosity")
    ax.set_ylabel("Size")
    ax.set_zlabel("Density")
    for i in range(len(materials)):
        ax.text(x[i], y[i], z[i], s[i])
    plt.show()


if __name__ == "__main__":
    main()
