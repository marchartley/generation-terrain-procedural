#include "Curve1D.h"

Curve1D::Curve1D()
{

}

Curve1D::Curve1D(const std::vector<Vector3>& points)
    : BSpline(points)
{
    std::sort(this->points.begin(), this->points.end(), [&](const Vector3& A, const Vector3& B) { return A.x < B.x; });
}

float Curve1D::get(float x) const
{
    for (size_t i = 0; i < this->size() - 1; i++) {
        auto& prev = this->points[i];
        auto& next = this->points[i + 1];
        if (prev.x <= x && x <= next.x) {
            float prevX = 0.f;
            float nextX = next.x - prev.x;
            x -= prev.x;
            if (prevX == nextX) return prev.y; // Strange case of 2 points with same X value
            float t = x / nextX;
            return prev.y * (1.f - t) + next.y * t; // Linear interpolation
        }
    }
    throw std::domain_error("[Curve1D error]: Trying to get value for x = " + std::to_string(x) + " but curve is defined on [" + std::to_string(points.front().x) + ", " + std::to_string(points.back().x) + "\nCurve is defined as " + this->toString());
}
