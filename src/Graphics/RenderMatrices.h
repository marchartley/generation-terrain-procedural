#ifndef RENDERMATRICES_H
#define RENDERMATRICES_H
#include "DataStructure/Matrix3.h"
#include <cmath>
#include <stdexcept>

// Matrix3 stores x contiguously. Here x is the row and y is the column, so
// data.data() has the column-major layout expected by GL_FALSE matrix uploads.
// Matrix3's ordinary operator* is elementwise; transform products are explicit.
namespace RenderMatrices {
template<class T>
Matrix3<T> multiply(const Matrix3<T>& a, const Matrix3<T>& b)
{
    Matrix3<T> result(4,4,1,T(0));
    for (int row=0; row<4; ++row)
        for (int col=0; col<4; ++col)
            for (int k=0; k<4; ++k)
                result(row,col) += a(row,k)*b(k,col);
    return result;
}

template<class T>
GridF floats(const Matrix3<T>& matrix)
{
    GridF result(4,4,1,0.f);
    for (size_t i=0;i<16;++i) result[i]=float(matrix[i]);
    return result;
}

inline GridF inverseTranspose(const float* columnMajor)
{
    Matrix3<double> a(4,4,1,0.0);
    auto inverse=Matrix3<double>::identity(4,4);
    for (int i=0;i<16;++i) a[i]=columnMajor[i];
    // Straightforward pivoted Gauss-Jordan elimination; optimization is deferred.
    for (int col=0;col<4;++col) {
        int pivot=col;
        for (int row=col+1;row<4;++row)
            if (std::abs(a(row,col))>std::abs(a(pivot,col))) pivot=row;
        if (std::abs(a(pivot,col))<1e-15) throw std::invalid_argument("Singular render matrix");
        for(int k=0;k<4;++k) {
            std::swap(a(col,k),a(pivot,k));
            std::swap(inverse(col,k),inverse(pivot,k));
        }
        const double divisor=a(col,col);
        for(int k=0;k<4;++k) { a(col,k)/=divisor; inverse(col,k)/=divisor; }
        for(int row=0;row<4;++row) {
            if(row==col)continue;
            const double factor=a(row,col);
            for(int k=0;k<4;++k) { a(row,k)-=factor*a(col,k); inverse(row,k)-=factor*inverse(col,k); }
        }
    }
    GridF result(4,4,1,0.f);
    for(int row=0;row<4;++row)
        for(int col=0;col<4;++col) result(row,col)=float(inverse(col,row));
    return result;
}

inline double length(const Vec3<double>& v) { return std::sqrt(v.dot(v)); }
inline Vec3<double> unit(const Vec3<double>& v) { return v/length(v); }

inline Matrix3<double> lookAt(const Vec3<double>& eye, const Vec3<double>& center, const Vec3<double>& up)
{
    const auto forward=unit(center-eye);
    const auto side=unit(forward.cross(up));
    const auto vertical=side.cross(forward);
    auto view=Matrix3<double>::identity(4,4);
    for(int col=0;col<3;++col) {
        view(0,col)=side.data()[col]; view(1,col)=vertical.data()[col]; view(2,col)=-forward.data()[col];
    }
    view(0,3)=-side.dot(eye); view(1,3)=-vertical.dot(eye); view(2,3)=forward.dot(eye);
    return view;
}

inline Matrix3<double> frustum(double left,double right,double bottom,double top,double nearZ,double farZ)
{
    Matrix3<double> projection(4,4,1,0.0);
    projection(0,0)=2*nearZ/(right-left); projection(1,1)=2*nearZ/(top-bottom);
    projection(0,2)=(right+left)/(right-left); projection(1,2)=(top+bottom)/(top-bottom);
    projection(2,2)=-(farZ+nearZ)/(farZ-nearZ); projection(2,3)=-2*farZ*nearZ/(farZ-nearZ);
    projection(3,2)=-1;
    return projection;
}
}
#endif
