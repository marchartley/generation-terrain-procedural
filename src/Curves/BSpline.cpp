#include "BSpline.h"

#include "Curves.h"

BSpline::BSpline() {}
BSpline::BSpline(const std::vector<Vector3> &points)
    : points(points)
{}
BSpline::BSpline(const Curve &curve)
    : points(curve.getPath())
{}

std::vector<Vector3> BSpline::getPath(int numberOfPoints) const
{
    if (numberOfPoints < 0)
        return this->points;

    std::vector<Vector3> points;
    for (int i = 0; i < numberOfPoints; i++) {
        float t = float(i) / float(numberOfPoints - 1);
        points.push_back(getPoint(t));
        if (!points.back().isValid()) {
            std::cout << "Error on t=" << t << std::endl;
        }
    }
    return points;
}

Vector3 BSpline::getPoint(float _x) const
{
    /*float res = 1 / (float)(numPoints() - 3);
    int iFloor = int(x / res);
    int iCeil = iFloor + 1;
    float resFloor = iFloor * res;
    float resCeil = iCeil * res;
    float x_prime = map(x, resFloor, resCeil, 0.f, 1.f);
    iFloor ++;
    iCeil ++;
    if (x >= 1.f) {
        iFloor --;
        iCeil --;
        x_prime = 1.f;
    }*/
    float x = std::clamp(_x, 0.f, 1.f);

    std::vector<Vector3> displayedPoints = this->points;
    size_t nbPoints = displayedPoints.size();

    if (nbPoints == 0) return Vector3::invalid;
    else if (nbPoints == 1) return displayedPoints[0];
    else if (nbPoints == 2) return displayedPoints[0] * (1.f - x) + displayedPoints[1] * x;

    float res = 1 / (float)(nbPoints - 1);
    int iFloor = int(x / res);
    int iCeil = int(x / res) + 1;
    float resFloor = iFloor * res;
    float resCeil = iCeil * res;
    float x_prime = map(x, resFloor, resCeil, 0.f, 1.f);

    if (!isClosed()) {
        displayedPoints.insert(displayedPoints.begin(), displayedPoints.front());
        displayedPoints.push_back(displayedPoints.back());
        if (_x < 1.f) {
            iFloor++;
            iCeil++;
        } else {
            x_prime = 1.f; // Small problem when t = 1.0
        }
        // nbPoints = displayedPoints.size();
    }

    auto P0 = displayedPoints[iFloor - 1];
    auto P1 = displayedPoints[iFloor];
    auto P2 = displayedPoints[iCeil];
    auto P3 = displayedPoints[iCeil + 1];

    return ((-P0 + P1 * 3.f - P2 * 3.f + P3) * x_prime * x_prime * x_prime
            + (P0 * 3.f - P1 * 6.f + P2 * 3.f) * x_prime * x_prime
            + (P0 * -3.f + P2 * 3.f) * x_prime
            + (P0 + P1 * 4.f + P2)) / 6.f;
}

Vector3 BSpline::getDerivative(float x, bool normalize) const
{
    return Vector3::invalid;
}

Vector3 BSpline::getSecondDerivative(float x, bool normalize) const
{
    return Vector3::invalid;
}

float BSpline::estimateClosestTime(const Vector3 &pos) const
{
    return toBezier(*this).estimateClosestTime(pos);
}

Vector3 BSpline::estimateClosestPos(const Vector3 &pos) const
{
    return getPoint(estimateClosestTime(pos));
}

float BSpline::estimateSqrDistanceFrom(const Vector3 &pos) const
{
    return (estimateClosestPos(pos) - pos).norm2();
}

float BSpline::length() const
{
    float length = 0;
    const auto& points = getPath(numPoints());
    for (size_t i = 0; i < points.size() - 1; i++) {
        length += (points[i] - points[i+1]).norm();
    }
    return length;
}

BSpline& BSpline::setPoint(int i, const Vector3 &newPos)
{
    this->points[pointIndex(i)] = newPos;
    return *this;
}

BSpline& BSpline::resamplePoints(int newNbPoints)
{
    if (newNbPoints < 0) newNbPoints = numPoints();
    std::vector<Vector3> newPoints(newNbPoints);
    for (int i = 0; i < newNbPoints; i++) {
        newPoints[i] = getPoint(float(i) / float(newNbPoints - 1));
    }
    this->points = newPoints;
    return *this;
}

BSpline& BSpline::reverseVertices()
{
    std::reverse(this->points.begin(), this->points.end());
    return *this;
}

std::pair<Vector3, Vector3> BSpline::AABBox() const
{
    Vector3 mini = Vector3::max();
    Vector3 maxi = Vector3::min();

    for (auto& p : points) {
        mini.x() = std::min(mini.x(), p.x());
        mini.y() = std::min(mini.y(), p.y());
        mini.z() = std::min(mini.z(), p.z());
        maxi.x() = std::max(maxi.x(), p.x());
        maxi.y() = std::max(maxi.y(), p.y());
        maxi.z() = std::max(maxi.z(), p.z());
    }
    return {mini, maxi};
}

BSpline& BSpline::scale(const Vector3 &factor)
{
    for (auto& p : points) {
        p *= factor;
    }
    return *this;
}

BSpline& BSpline::translate(const Vector3 &translation)
{
    for (auto& p : points) {
        p += translation;
    }
    return *this;
}

BSpline& BSpline::removeDuplicates()
{
    return *this;
}

size_t BSpline::numPoints() const
{
    return points.size();
}

std::string BSpline::toString() const
{
    std::ostringstream out;
    out << "B-Spline with " << this->points.size() << " points (" << (closed ? "closed" : "not closed") << ") :\n";
    for (auto& p : this->points)
        out << "- " << p << "\n";
    return out.str();
}

BSpline& BSpline::close()
{
    Curve::close();
    if (this->points.size() > 0 && this->points.front() != this->points.back()) {
        this->points.push_back(points[0]);
        this->points.push_back(points[1]);
        this->points.push_back(points[2]);
    }
    return *this;
}

Vector3& BSpline::operator[](size_t i)
{
    return this->points[(i + size()) % size()];
}
const Vector3& BSpline::operator[](size_t i) const
{
    return this->points[(i + size()) % size()];
}

