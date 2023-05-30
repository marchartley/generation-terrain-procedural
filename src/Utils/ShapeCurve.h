#ifndef SHAPECURVE_H
#define SHAPECURVE_H

#include "BSpline.h"

class ShapeCurve : public BSpline
{
public:
    ShapeCurve();
    ShapeCurve(std::vector<Vector3> points);
    ShapeCurve(BSpline path);

    bool contains(Vector3 pos, bool useNativeShape = true);
    float estimateDistanceFrom(Vector3 pos);
    float computeArea();
    float computeSignedArea();

    ShapeCurve& reverseVertices();

    ShapeCurve intersect(ShapeCurve other);

    Vector3 planeNormal();
    ShapeCurve grow(float increase);
    ShapeCurve shrink(float decrease);
    ShapeCurve& translate(Vector3 translation);

    ShapeCurve& removeDuplicates();

    std::vector<Vector3> closedPath();

    std::vector<Vector3> randomPointsInside(int numberOfPoints = 1);

    ShapeCurve& scale(float factor);

    ShapeCurve merge(ShapeCurve other);
};

struct ClipVertex {
    Vector3 coord;
    bool isEntry;
    bool isExit;
    bool isInside;
    ClipVertex* neighbor = nullptr;
    ClipVertex* prev = nullptr;
    ClipVertex* next = nullptr;
    int shapeID = -1;
    int index = -1;

    ClipVertex(Vector3 p) : ClipVertex(p, false, false, false) {}
    ClipVertex(Vector3 p, bool entry, bool exit, bool inside)
        : coord(p), isEntry(entry), isExit(exit), isInside(inside)
    {}
};
int getIndex(int index, size_t size);
int markEntriesExits(std::vector<ClipVertex*>& poly, bool currentlyInside, int shapeID);

#endif // SHAPECURVE_H
