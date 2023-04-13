import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import os.path
from typing import Tuple, List, Callable, Union, Optional, Any
from matplotlib.widgets import Button
from Vectors import Vector2D, Vector3D, line_intersection
import Vectors
from noise import perlin
from sketch_maker import LineBuilder1D, LineBuilder2D, LineBuilderRadial, SketchManagement, LineBuilder


def distanceMapFromSketches(islandSketch: SketchManagement) -> np.ndarray:
    dims = (100, 100, len(islandSketch.lineBuilders) + 1) # Add the "borders" as last dimension
    heightmap = np.zeros(dims)

    for _y in range(dims[0]):
        for _x in range(dims[1]):
            x = (_x / dims[1]) * 2.0 - 1.0  # Transferring on [-1, 1]^2
            y = (_y / dims[0]) * 2.0 - 1.0
            pos = Vector2D(x, y)
            for iChannel, line in enumerate(islandSketch.lineBuilders):
                minDist = max(dims)
                curve = line.getCurve()
                for i in range(len(curve) - 1):
                    minDist = min(minDist, Vectors.distanceToLine(pos, curve[i], curve[i + 1]))
                heightmap[dims[0] - _y - 1, _x, iChannel] = minDist
            distToBorder = min([abs(x), abs(y), abs(1-x), abs(1-y)])
            heightmap[dims[0] - _y - 1, _x, -1] = distToBorder
    return heightmap

def genAndSaveDistanceMap(islandSketch: SketchManagement):
    image = distanceMapFromSketches(islandSketch)
    image = image[:, :, :3]
    # image = 1 - np.clip(image * 5.0, 0, 1)
    plt.imsave("test.png", image)

def splitProfileOnMarkers(profileSketch: LineBuilder, islandSketches: SketchManagement, sliceCut: Tuple[Vector2D, Vector2D]) -> List[Tuple[int, int, List[float]]]:
    """Extract the curves made by the profile depending on the distance between each island sketch"""
    markers: List[Tuple[float, int]] = []

    sliceDirection = (sliceCut[1] - sliceCut[0]).normalize()
    for sketchID, sketch in enumerate(islandSketches.lineBuilders):
        intersections = sketch.intersection(*sliceCut)
        for intersect in intersections:
            # distanceOnSlice = sliceDirection.dot(intersect - sliceCut[0])
            # print()
            distanceOnSlice = intersect.x
            markers.append((distanceOnSlice, sketchID))

    markers.sort(key = lambda a: a[0])

    borderMarker = len(islandSketches.lineBuilders)

    profile = profileSketch.getCurve()
    curves: List[Tuple[int, int, List[float]]] = []
    begin = (-1.0, borderMarker)  # markers.pop(0)
    ending = markers[0]

    startingIndex = 0
    for i, p in enumerate(profile):
        x, y = p.x, p.y
        if x >= ending[0]:  # End of the curve, beginning of a new one
            curves.append((begin[1], ending[1], [_p.y for _p in profile[startingIndex:i]]))
            if len(markers) > 0:
                startingIndex = i
                begin = markers.pop(0)
                if len(markers) == 0:
                    ending = (1.0, borderMarker)
                else:
                    ending = markers[0]
    return curves

def interpolateOnCurve(curve: List[Any], t: float) -> Vector2D:
    if t >= 1.0:
        return curve[-1]
    i = t * len(curve)
    a = i - math.floor(i)
    if math.floor(i) < 0 or math.ceil(i) >= len(curve):
        b = 0
    return curve[math.floor(i)] * a + curve[math.ceil(i)] * (1 - a)

def heightmapFromSketches(profileSketch: LineBuilder, islandSketches: SketchManagement, sliceCut: Tuple[Vector2D, Vector2D]):
    sequences = splitProfileOnMarkers(profileSketch, islandSketches, sliceCut)
    distances = distanceMapFromSketches(islandSketches)
    heights = np.zeros((distances.shape[0], distances.shape[1]))

    for _y in range(heights.shape[0]):
        for _x in range(heights.shape[1]):
            closestFeatures = [i[0] for i in sorted(enumerate(distances[_y, _x, :-1]), key=lambda x:x[1])]  # Get indices of sorted array
            top2 = closestFeatures[:2]
            valid_sequences = [s[2] for s in sequences if list(s[:2]) == list(top2)]
            if valid_sequences:
                sequence = valid_sequences[0]
                dist0 = distances[_y, _x, top2[0]]
                dist1 = distances[_y, _x, top2[1]]
                if len(sequence) > 1 and dist0 + dist1 != 0.0:
                    heights[_y, _x] = interpolateOnCurve(sequence, dist0 / (dist0 + dist1))
                else:
                    heights[_y, _x] = sequence[0]
    return heights

def genAndSaveHeightMap(profileSketch: LineBuilder, islandSketches: SketchManagement, sliceCut: Tuple[Vector2D, Vector2D]):
    path = "test.png"
    image = heightmapFromSketches(profileSketch, islandSketches, sliceCut)
    # image /= np.max(image)
    rgb = np.zeros((image.shape[0], image.shape[1], 3))
    rgb[:, :, 0] = rgb[:, :, 1] = rgb[:, :, 2] = image
    plt.imsave(path, rgb)
    print("Heightmap saved at " + os.path.abspath(path))

def main():
    sketch_names = ["Reef", "Island", "Passes"]
    sketch_colors = ["orange", "green", "blue"]
    sketch_types = [LineBuilderRadial, LineBuilderRadial, LineBuilder2D]
    waterLevel = 0.7

    # Creates a top view for island sketching and a side view for profile editing
    fig, axes = plt.subplots(1, 2, squeeze=False)
    topViewAx: plt.Axes = axes[0, 0]
    sideViewAx: plt.Axes = axes[0, 1]
    fig.subplots_adjust(bottom=0.2)

    topViewAx.set_title('Top view')
    topViewAx.set_xlim(-1, 1)
    topViewAx.set_ylim(-1, 1)
    sideViewAx.set_title('Side view')
    sideViewAx.set_xlim(-1, 1)
    sideViewAx.set_ylim(0, 1)

    # Represent the center of the map with an ellipse
    # topViewAx.add_patch(matplotlib.patches.Circle((0, 0), 0.1))

    # Add sketching for top island sketching
    islandSketches = SketchManagement(topViewAx)

    buttons: List[Button] = []
    button_size = 0.5 / len(sketch_names)
    for i_sketch in range(len(sketch_names)):
        islandSketches.addSketch(color = sketch_colors[i_sketch], sketch_type = sketch_types[i_sketch])
        sketchID = i_sketch
        ax_button = fig.add_axes([0.5 + button_size * sketchID, 0.05, button_size - 0.01, 0.075])
        buttons.append(Button(ax_button, sketch_names[sketchID]))
        def activationFunction(id):
            def action(event):
                islandSketches.activate(id)
            return action
        buttons[-1].on_clicked(activationFunction(sketchID))

    # Add button for generating a distance map / heightmap
    ax_button = fig.add_axes([0.0, 0.05, 0.09, 0.075])
    distance_button = Button(ax_button, "Gen height map")
    distance_button.on_clicked(lambda e: genAndSaveHeightMap(profileSketching.lineBuilders[0], islandSketches, sliceCut))
    # Add button for splitting sequences
    ax_button = fig.add_axes([0.1, 0.05, 0.09, 0.075])
    splitting_button = Button(ax_button, "Split profile")
    splitting_button.on_clicked(lambda e: splitProfileOnMarkers(profileSketching.lineBuilders[0], islandSketches, sliceCut))

    # Profile will be defined by the horizontal plane passing through center
    sliceCut = (Vector2D(-1, 0), Vector2D(1, 0))
    topViewAx.plot([sliceCut[0].x, sliceCut[1].x], [sliceCut[0].y, sliceCut[1].y], linestyle="--", color="blue")

    # Each time the top view sketch is modified, reposition markers on the side view
    def updateSideViewMarkers():
        for line in sideViewAx.lines:
            line.set_data([], [])
        # Add water level to profile editing
        sideViewAx.axhline(y = waterLevel, color = "blue")
        for sketch in islandSketches.lineBuilders:
            intersections = sketch.intersection(*sliceCut)
            for intersect in intersections:
                sideViewAx.axvline(x = intersect.x, color = sketch.color)
        islandSketches.draw()
        profileSketching.draw()
        fig.canvas.draw()
        # fig.canvas.flush_events()

    # Add sketching for profile editing
    profileSketching: SketchManagement = SketchManagement(sideViewAx)
    profileSketching.addSketch(line = LineBuilder1D(sideViewAx, 20, "blue"))

    # Testing part
    def centeredCircle(radius:float, randomness: float = 0.2) -> List[Vector2D]:
        noise = perlin.SimplexNoise(10)
        points: List[Vector2D] = []
        nbPoints = 30
        for i in range(nbPoints):
            angle = i * math.tau / (nbPoints - 1)
            vertexUnitPos = Vector2D(math.cos(angle), math.sin(angle))
            noiseValue = noise.noise2(vertexUnitPos.x, vertexUnitPos.y) * randomness
            points.append(vertexUnitPos * (radius * (1 - noiseValue)))
        return points
    randomCoralCurve = centeredCircle(0.5)
    islandSketches.lineBuilders[0].setCurve(randomCoralCurve)
    randomIslandCurve = centeredCircle(0.25)
    islandSketches.lineBuilders[1].setCurve(randomIslandCurve)

    _randomProfileCurve = [.0, .5, .6, .7, .7, .6, .6, .9, .9, .6, .6, .7, .7, .6, .5, .0]
    randomProfileCurve = []
    for i in range(len(_randomProfileCurve)):
        randomProfileCurve.append(Vector2D((2 * i/(len(_randomProfileCurve) - 1)) - 1, _randomProfileCurve[i]))
    profileSketching.lineBuilders[0].setCurve(randomProfileCurve)

    islandSketches.onChange(updateSideViewMarkers)
    # Todo:
    #  - From profile, get heightmap depending on coral position and island position
    #  - Add passes in the process
    #  - Create heightmap from input and distance map
    #  - Profile differs depending on "r" and "d" (island width, island-coral distance)

    updateSideViewMarkers()
    genAndSaveHeightMap(profileSketching.lineBuilders[0], islandSketches, sliceCut)

    plt.show()


if __name__ == "__main__":
    main()
