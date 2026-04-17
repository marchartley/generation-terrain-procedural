#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>

class Matrix // : public std::vector<std::vector<float>>
{
public:
    Matrix();
    Matrix(size_t n, size_t m, float default_value = 0.f);
    Matrix(size_t n, size_t m, const float* data);
    // Matrix(size_t n, size_t m, const float** data); // Never used...
    Matrix(size_t n, size_t m, const std::vector<float>& data);
    Matrix(const std::vector<std::vector<float>>& data);

    float det() const;
    Matrix adj() const;
    Matrix inverse() const;
    Matrix cofactors() const;
    Matrix submatrix(size_t rowToIgnore, size_t colToIgnore) const;
    Matrix transpose() const;
    Matrix product(const Matrix& m) const;
    float trace() const;

    const size_t& rows() const { return sizeY; } //this->size(); }
    const size_t& cols() const { return sizeX; } // (*this)[0].size(); }

    Matrix col(size_t colIndex) const;
    Matrix row(size_t rowIndex) const;

    Matrix abs() const;

    float maxCoeff() const;

    Matrix leftCols(size_t nbCols) const;
    Matrix rightCols(size_t nbCols) const;

    static Matrix identity(size_t size);

    std::string displayValues() const;
    std::string toString() const;
    std::string displayValuesOneLine() const;

    inline float& operator()(size_t x, size_t y) { return this->unsafe(x, y); }
    inline const float& operator()(size_t x, size_t y) const { return this->unsafe(x, y); }
    inline float& unsafe(size_t x, size_t y) { return _data[y * sizeX + x]; } //return (*this)[x][y]; }
    inline const float& unsafe(size_t x, size_t y) const { return _data[y * sizeX + x]; } // return (*this)[x][y]; }
    inline const std::vector<float>& data() const { return _data; }

    size_t size() const { return data().size(); }

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

    static std::pair<Matrix, Matrix> gramSchmidtQR(const Matrix& A);
    static std::vector<float> backSubstitution(const Matrix& R, const std::vector<float>& b);
    static std::vector<float> solve(const Matrix& A, const Matrix& b);
    static std::vector<float> solve(const Matrix& A, const std::vector<float>& b);


    size_t sizeX;
    size_t sizeY;
    std::vector<float> _data;
};

#endif // MATRIX_H
