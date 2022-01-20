#include "Utils/BSpline.h"

BSpline::BSpline()
{

}
BSpline::BSpline(std::vector<Vector3> points) : points(points)
{

}
BSpline::BSpline(int numberOfPoints) {
    for(int i = 0; i < numberOfPoints; i++) {
        this->points.push_back(Vector3::random());
    }
}

std::vector<Vector3> BSpline::getPath(float resolution)
{
    std::vector<Vector3> path;
    for (float x = 0.0; x < 1.0; x += resolution)
        path.push_back(this->getPoint(x));
    path.push_back(this->getPoint(1.0));
    return path;
}

Vector3 BSpline::getPoint(float x)
{
    if (points.size() >= 4)
        return this->getCatmullPoint(x);

    if(points.size() == 0)
        return Vector3();

    std::vector<Vector3> controls = points;
    while (controls.size() > 1)
    {
        std::vector<Vector3> newCtrls;
        for (size_t i = 0; i < controls.size() - 1; i++) {
            newCtrls.push_back(this->getPoint(x, controls[i], controls[i+1]));
        }
        controls = newCtrls;
    }
    return controls[0];
}
Vector3 BSpline::getPoint(float x, Vector3 a, Vector3 b)
{
    return a * (1 - x) + b * x;
}

float CatmullNextT(Vector3 P0, Vector3 P1, float t_prev, float alpha)
{
    return std::pow((P0 - P1).norm2(), alpha) + t_prev;
}
template <class T>
T map(T x, T prev_min, T prev_max, T new_min, T new_max)
{
    return ((x - prev_min) / (prev_max - prev_min)) * (new_max - new_min) + new_min;
}

Vector3 BSpline::getCatmullPoint(float x)
{
    float alpha = 0.5f;

    //alpha /= 2.f;

    float res = 1 / (float)(this->points.size() - 1);
    int iFloor = int(x / res);
    int iCeil = int(x / res) + 1;
    float resFloor = iFloor * res;
    float resCeil = iCeil * res;
    float x_prime = map(x, resFloor, resCeil, 0.f, 1.f);

    if (iFloor <= 0) {
        return this->getPoint(x_prime, this->points[0], this->points[1]);
    } else if (iCeil >= int(this->points.size()) - 1) {
        return this->getPoint(x_prime, this->points[iFloor], this->points[iCeil]);
    }

    Vector3 P0 = points[iFloor - 1];
    Vector3 P1 = points[iFloor - 0];
    Vector3 P2 = points[iCeil + 0];
    Vector3 P3 = points[iCeil + 1];

    float t0 = 0;
    float t1 = CatmullNextT(P0, P1, t0, alpha);
    float t2 = CatmullNextT(P1, P2, t1, alpha);
    float t3 = CatmullNextT(P2, P3, t2, alpha);

    float t = map(x_prime, 0.f, 1.f, t1, t2);

    Vector3 A1 = P0 * (t1 - t) / (t1 - t0) + P1 * (t - t0) / (t1 - t0);
    Vector3 A2 = P1 * (t2 - t) / (t2 - t1) + P2 * (t - t1) / (t2 - t1);
    Vector3 A3 = P2 * (t3 - t) / (t3 - t2) + P3 * (t - t2) / (t3 - t2);

    Vector3 B1 = A1 * (t2 - t) / (t2 - t0) + A2 * (t - t0) / (t2 - t0);
    Vector3 B2 = A2 * (t3 - t) / (t3 - t1) + A3 * (t - t1) / (t3 - t1);

    Vector3 C  = B1 * (t2 - t) / (t2 - t1) + B2 * (t - t1) / (t2 - t1);
    return C;
}
