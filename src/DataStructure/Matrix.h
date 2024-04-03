#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>

class Matrix : public std::vector<std::vector<float>>
{
public:
    Matrix();
    Matrix(int n, int m, float default_value = 0.f);
    Matrix(int n, int m, float* data);
    Matrix(int n, int m, float** data);
    Matrix(const std::vector<std::vector<float>>& data);

    float det() const;
    Matrix adj() const;
    Matrix inverse() const;
    Matrix cofactors() const;
    Matrix submatrix(size_t rowToIgnore, size_t colToIgnore) const;
    Matrix transpose() const;
    Matrix product(const Matrix& m) const;
    float trace() const;

    int rows() const { return this->size(); }
    int cols() const { return (*this)[0].size(); }

    static Matrix identity(int size);

    std::string displayValues() const;
    std::string toString() const;
    std::string displayValuesOneLine() const;

    float* toArray(bool rowByRow = true) const;
    std::vector<float> toStdVector(bool rowByRow = true) const;


//    friend Matrix operator+(Matrix a, Matrix& b);
    friend Matrix operator+(const Matrix& a, const Matrix& b);
    Matrix& operator+=(const Matrix& o);
    friend Matrix operator-(const Matrix& a, const Matrix& b);
    Matrix& operator-=(const Matrix& o);
    friend Matrix operator*(const Matrix& a, const Matrix& o);
    Matrix& operator*=(const Matrix& o);
    Matrix operator/(const Matrix& o);
    Matrix& operator/=(const Matrix& o);
    Matrix operator*(float o);
    Matrix& operator*=(float o);
    Matrix operator/(float o);
    Matrix& operator/=(float o);
    Matrix operator+(float o);
    Matrix& operator+=(float o);
    Matrix operator-(float o);
    Matrix& operator-=(float o);

    static Matrix matprod(const Matrix& A, const Matrix& B);

    Matrix toHomogeneous() const;

    friend std::ostream& operator<<(std::ostream& io, const Matrix& m);
    friend std::ostream& operator<<(std::ostream& io, Matrix* m);
};

#endif // MATRIX_H
