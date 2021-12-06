#include "BSpline.h"

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

Vector3 BSpline::getPoint(float x)
{
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
