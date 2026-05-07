#ifndef SHAPECURVE_H
#define SHAPECURVE_H

#include "Curves.h"

class Contour // : public CatmullRomSpline
{
public:
    Contour();
    Contour(std::vector<Vector3> points);
    Contour(const Curve& path);
    Contour(Curve* path);
    Contour(std::shared_ptr<Curve> path);
    ~Contour();

    bool contains(const Vector3& pos, bool useNativeShape = true) const;
    bool containsXY(const Vector3& pos, bool useNativeShape = true, int increaseAccuracy = 0) const;
    float estimateDistanceFrom(const Vector3& pos) const;
    float estimateSignedDistanceFrom(const Vector3& pos) const;
    float computeArea() const;
    float computeSignedArea() const;

    Vector3 centroid() const;

    // Contour intersect(Contour other);

    Vector3 planeNormal() const;
    Contour grow(float increase);
    Contour shrink(float decrease);
    Contour& translate(const Vector3& translation);

    Contour& removeDuplicates();

    std::vector<Vector3> closedPath() const;

    std::vector<Vector3> randomPointsInside(int numberOfPoints = 1);

    Contour merge(Contour other);

    Contour& resamplePoints(int newNbPoints = -1);

    Contour& setPoint(int i, const Vector3 &newPos);

    Polyline computeConvexHull() const;

    static std::vector<Vector3> circle(float radius, const Vector3& center, int nbPoints);

public:
    std::shared_ptr<Curve> curve;
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

    ClipVertex(const Vector3& p) : ClipVertex(p, false, false, false) {}
    ClipVertex(const Vector3& p, bool entry, bool exit, bool inside)
        : coord(p), isEntry(entry), isExit(exit), isInside(inside)
    {}
};
int getIndex(int index, size_t size);
int markEntriesExits(std::vector<ClipVertex*>& poly, bool currentlyInside, int shapeID);



std::vector<float> computeGreenCoordinates(const Vector3& p, const Contour& polygon);
Vector3 computePointFromGreenCoordinates(const std::vector<float>& greenCoords, const Contour& polygon);

Vector3 randomPointInTriangle(const Vector3& A, const Vector3& B, const Vector3& C);

#endif // SHAPECURVE_H
