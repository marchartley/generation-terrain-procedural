#ifndef PLOTTINGUTILS_H
#define PLOTTINGUTILS_H

#include "DataStructure/Matrix3.h"

#include "Utils/Signals.h"

struct PlottingUtils {
    template <class T>
    static Matrix3<T>& drawLine(Matrix3<T>& img, const T& color, const Vector3& start, const Vector3& end, int strokeWidth = 1);
    template <class T>
    static Matrix3<T>& drawCircle(Matrix3<T>& img, const T& color, const Vector3& center, float radius);
    template <class T>
    static Matrix3<T>& drawCurve(Matrix3<T>& img, const T& color, const CatmullRomSpline& spline, int strokeWidth = 1);
};

template <class T>
Matrix3<T>& PlottingUtils::drawLine(Matrix3<T>& img, const T& color, const Vector3& start, const Vector3& end, int strokeWidth) {
    auto line = (end - start);
    int dx = line.x();
    int dy = line.y();

    // calculate steps required for generating pixels
    int steps = abs(dx) > abs(dy) ? abs(dx) : abs(dy);

    // calculate increment in x & y for each steps
    float Xinc = dx / (float)steps;
    float Yinc = dy / (float)steps;

    // Put pixel for each step
    auto p = start;
    for (int i = 0; i <= steps; i++) {
        Vector3i _p = p;
        for (int off = -strokeWidth / 2; off <= strokeWidth / 2; off++) {
            img[_p + Vector3i(off, 0, 0)] = color;
            img[_p + Vector3i(0, off, 0)] = color;
        }
        p.x() += Xinc;
        p.y() += Yinc;
    }

    return img;
}


template<class T>
Matrix3<T>& PlottingUtils::drawCircle(Matrix3<T> &img, const T &color, const Vector3 &center, float radius)
{
    auto drawCircle = [&](int xc, int yc, int x, int y){
        img(Vector3i(xc + x, yc + y)) = color;
        img(Vector3i(xc - x, yc + y)) = color;
        img(Vector3i(xc + x, yc - y)) = color;
        img(Vector3i(xc - x, yc - y)) = color;
        img(Vector3i(xc + y, yc + x)) = color;
        img(Vector3i(xc - y, yc + x)) = color;
        img(Vector3i(xc + y, yc - x)) = color;
        img(Vector3i(xc - y, yc - x)) = color;
    };

    int xc = center.x();
    int yc = center.y();
    int r = radius;

    // using Bresenham's algorithm
    int x = 0, y = r;
    int d = 3 - 2 * r;
    drawCircle(xc, yc, x, y);
    while (y >= x){

        if (d > 0) {
            y--;
            d = d + 4 * (x - y) + 10;
        }
        else
            d = d + 4 * x + 6;
        x++;
        drawCircle(xc, yc, x, y);
    }
    return img;
}


template <class T>
Matrix3<T>& PlottingUtils::drawCurve(Matrix3<T>& img, const T& color, const CatmullRomSpline& spline, int strokeWidth)
{
    int res = 100;
    for (int i = 0; i < res + 1; i++) {
        float t0 = float(i) / float(res + 1);
        float t1 = float((i < res + 1 ? i + 1 : 0)) / float(res + 1);
        drawLine(img, color, spline.getPoint(t0), spline.getPoint(t1), strokeWidth);
    }
    return img;
}



#endif // PLOTTINGUTILS_H
