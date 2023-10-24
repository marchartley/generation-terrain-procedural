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
    Matrix(std::vector<std::vector<float>> data);

    float det() const;
    Matrix adj() const;
    Matrix inverse() const;
    Matrix cofactors() const;
    Matrix submatrix(size_t rowToIgnore, size_t colToIgnore) const;
    Matrix transpose() const;
    Matrix product(Matrix m) const;
    float trace() const;

    std::string displayValues() const;
    std::string toString() const;
    std::string displayValuesOneLine() const;


//    friend Matrix operator+(Matrix a, Matrix& b);
    friend Matrix operator+(Matrix a, Matrix b);
    Matrix& operator+=(const Matrix& o);
    friend Matrix operator-(Matrix a, const Matrix& b);
    Matrix& operator-=(const Matrix& o);
    friend Matrix operator*(Matrix a, const Matrix& o);
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

    static Matrix matprod(Matrix A, Matrix B);

    friend std::ostream& operator<<(std::ostream& io, const Matrix& m);
    friend std::ostream& operator<<(std::ostream& io, Matrix* m);
};

#endif // MATRIX_H
