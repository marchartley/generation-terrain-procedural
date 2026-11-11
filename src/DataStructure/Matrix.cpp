#include "DataStructure/Matrix.h"


#include <cmath>
#include <sstream>
#include <limits>
#include <algorithm>
#include <stdexcept>

namespace {
size_t count(size_t x, size_t y) {
    if (y && x > std::vector<float>().max_size() / y)
        throw std::length_error("Matrix size overflow");
    return x * y;
}
void validate(const Matrix &m) {
    if (count(m.sizeX, m.sizeY) != m.size())
        throw std::invalid_argument("Matrix payload does not match dimensions");
}
void square(const Matrix &m) {
    validate(m);
    if (m.cols() != m.rows())
        throw std::invalid_argument("Matrix must be square");
}
} // namespace

Matrix::Matrix()
    : Matrix(1, 1, 0.0)
{}

Matrix::Matrix(size_t n, size_t m, float default_value)
{
    sizeX = n;
    sizeY = m;
    _data.assign(count(n, m), default_value);
}
Matrix::Matrix(size_t n, size_t m, const float* data){
    sizeX = n;
    sizeY = m;
    auto size = count(n, m);
    if (size && !data)
        throw std::invalid_argument("Null matrix payload");
    if (size)
        _data.assign(data, data + size);
}

Matrix::Matrix(size_t n, size_t m, const std::vector<float> &data)
    : Matrix(n, m)
{
    if (data.size() != count(n, m))
        throw std::invalid_argument("Matrix payload size mismatch");
    _data = data;
}
Matrix::Matrix(const std::vector<std::vector<float> > &data)
{
    sizeX = data.empty() ? 0 : data.front().size();
    sizeY = data.size();
    _data.reserve(count(sizeX, sizeY));
    for (const auto &row : data) {
        if (row.size() != sizeX)
            throw std::invalid_argument("Ragged matrix");
        _data.insert(_data.end(), row.begin(), row.end());
    }
}

float Matrix::det() const
{
    square(*this);
    const size_t n = rows();
    std::vector<double> a(_data.begin(), _data.end());
    double determinant = 1;
    for (size_t k = 0; k < n; ++k) {
        size_t pivot = k;
        for (size_t r = k + 1; r < n; ++r)
            if (std::abs(a[r * n + k]) > std::abs(a[pivot * n + k]))
                pivot = r;
        if (a[pivot * n + k] == 0)
            return 0;
        if (pivot != k) {
            for (size_t c = 0; c < n; ++c)
                std::swap(a[k * n + c], a[pivot * n + c]);
            determinant = -determinant;
        }
        double d = a[k * n + k];
        determinant *= d;
        for (size_t r = k + 1; r < n; ++r) {
            double f = a[r * n + k] / d;
            for (size_t c = k + 1; c < n; ++c)
                a[r * n + c] -= f * a[k * n + c];
        }
    }
    return float(determinant);
}
Matrix Matrix::adj() const
{
    return this->cofactors().transpose();
}
Matrix Matrix::inverse() const
{
    square(*this);
    size_t n = rows();
    Matrix result = identity(n);
    std::vector<double> a(_data.begin(), _data.end()), inv(result._data.begin(), result._data.end());
    for (size_t k = 0; k < n; ++k) {
        size_t pivot = k;
        for (size_t r = k + 1; r < n; ++r)
            if (std::abs(a[r * n + k]) > std::abs(a[pivot * n + k]))
                pivot = r;
        double d = a[pivot * n + k];
        if (d == 0 || !std::isfinite(d))
            throw std::domain_error("Singular or nonfinite matrix");
        if (pivot != k)
            for (size_t c = 0; c < n; ++c) {
                std::swap(a[k * n + c], a[pivot * n + c]);
                std::swap(inv[k * n + c], inv[pivot * n + c]);
            }
        for (size_t c = 0; c < n; ++c) {
            a[k * n + c] /= d;
            inv[k * n + c] /= d;
        }
        for (size_t r = 0; r < n; ++r)
            if (r != k) {
                double f = a[r * n + k];
                for (size_t c = 0; c < n; ++c) {
                    a[r * n + c] -= f * a[k * n + c];
                    inv[r * n + c] -= f * inv[k * n + c];
                }
            }
    }
    for (double value : inv)
        if (!std::isfinite(value) || std::abs(value) > std::numeric_limits<float>::max())
            throw std::overflow_error("Matrix inverse cannot be represented as float");
    std::copy(inv.begin(), inv.end(), result._data.begin());
    return result;
}
Matrix Matrix::cofactors() const
{
    square(*this);
    Matrix result(cols(), rows());
    for (size_t r = 0; r < rows(); ++r)
        for (size_t c = 0; c < cols(); ++c)
            result.unsafe(c, r) = ((r + c) % 2 ? -1.f : 1.f) * submatrix(r, c).det();
    return result;
}
Matrix Matrix::submatrix(size_t rowToIgnore, size_t colToIgnore) const
{
    validate(*this);
    if (rowToIgnore >= rows() || colToIgnore >= cols())
        throw std::out_of_range("Invalid submatrix exclusion");
    Matrix result(cols() - 1, rows() - 1);
    size_t i = 0;
    for (size_t r = 0; r < rows(); ++r)
        if (r != rowToIgnore)
            for (size_t c = 0; c < cols(); ++c)
                if (c != colToIgnore)
                    result._data[i++] = unsafe(c, r);
    return result;
}

Matrix Matrix::transpose() const
{
    Matrix temp(rows(), cols());
    for(size_t row = 0; row < rows(); row++)
        for(size_t col = 0; col < cols(); col++)
            temp.unsafe(row, col) = this->unsafe(col, row);
    return temp;
}

Matrix Matrix::product(const Matrix& mat) const
{
    validate(*this);
    validate(mat);
    size_t m = this->cols();
    size_t l = this->rows();
    size_t n = mat.cols();
    size_t m2 = mat.rows();
    if (m != m2) {
        std::ostringstream oss;
        oss << "Error on product between matrices. A defined as " << m << "x" << l << " and B defined as " << n << "x" << m2; // << ":\n";
        // oss << this->toString() << "*\n" << mat.toString();
        throw std::domain_error(oss.str());
    }
    Matrix temp(n, l);

    for (size_t row = 0; row < l; ++row) {
        for (size_t k = 0; k < m; ++k) {
            const auto a = this->unsafe(k, row);
            for (size_t col = 0; col < n; ++col) {
                temp.unsafe(col, row) += a * mat.unsafe(col, k);
            }
        }
    }
    return temp;
}

float Matrix::trace() const
{
    if (sizeX != sizeY)
        throw std::domain_error("Error on matrix trace computation. Matrix is not square (" + std::to_string(sizeX) + " x " + std::to_string(sizeY) + ")");
    float trace = 0.0;
    for (size_t i = 0; i < cols(); i++)
        trace += this->unsafe(i, i);
    return trace;
}

Matrix Matrix::col(size_t colIndex) const
{
    validate(*this);
    if (colIndex >= cols())
        throw std::out_of_range("Column index");
    Matrix res(1, this->rows());
    for (size_t i = 0; i < this->rows(); i++) {
        res._data[i] = this->unsafe(colIndex, i);
    }
    return res;
}

Matrix Matrix::row(size_t rowIndex) const
{
    validate(*this);
    if (rowIndex >= rows())
        throw std::out_of_range("Row index");
    Matrix res(this->cols(), 1);
    for (size_t i = 0; i < this->cols(); i++) {
        res._data[i] = this->unsafe(i, rowIndex);
    }
    return res;
}

Matrix Matrix::abs() const
{
    Matrix res = *this;
    for (size_t i = 0; i < res._data.size(); i++) {
        res._data[i] = std::abs(res._data[i]);
    }
    return res;
}

float Matrix::maxCoeff() const
{
    float max = std::numeric_limits<float>::lowest();

    for (size_t i = 0; i < this->size(); i++)
        max = std::max(max, this->_data[i]);
    return max;
}

Matrix Matrix::leftCols(size_t nbCols) const
{
    validate(*this);
    if (nbCols > cols())
        throw std::out_of_range("Invalid column count");
    Matrix result(nbCols, rows());
    for (size_t r = 0; r < rows(); ++r)
        for (size_t c = 0; c < nbCols; ++c)
            result.unsafe(c, r) = unsafe(c, r);
    return result;
}

Matrix Matrix::rightCols(size_t nbCols) const
{
    validate(*this);
    if (nbCols > cols())
        throw std::out_of_range("Invalid column count");
    Matrix result(nbCols, rows());
    for (size_t r = 0; r < rows(); ++r)
        for (size_t c = 0; c < nbCols; ++c)
            result.unsafe(c, r) = unsafe(cols() - nbCols + c, r);
    return result;
}

Matrix Matrix::identity(size_t size)
{
    Matrix m(size, size);
    for (size_t i = 0; i < size; i++) {
        m(i, i) = 1.f;
    }
    return m;
}

std::string Matrix::displayValues() const
{
    std::string txt = "";
    for (size_t row = 0; row < rows(); row ++) {
        for (size_t col = 0; col < cols(); col ++) {
            txt += std::to_string(this->unsafe(col, row)) + "\t";
        }
        txt += "\n";
    }
    return txt;
}

std::string Matrix::toString() const
{
    return "Matrix (" + std::to_string(cols()) + "x" + std::to_string(rows()) + ") :\n" + this->displayValues();
    //    return txt;
}

std::string Matrix::displayValuesOneLine() const
{
    std::string txt = "";
    for (size_t row = 0; row < rows(); row ++) {
        for (size_t col = 0; col < cols(); col ++) {
            txt += std::to_string(this->unsafe(col, row)) + " ";
        }
    }
    return txt;
}

Matrix operator+(const Matrix& a, const Matrix& b)
{
    Matrix temp = a;
    temp += b;
    return temp;
}
Matrix& Matrix::operator+=(const Matrix& o)
{
    if (sizeX != o.sizeX || sizeY != o.sizeY)
        throw std::domain_error("Error on matrix addition. Matrix A (" + std::to_string(sizeX) + "x" + std::to_string(sizeY) + ") does not match Matrix B (" + std::to_string(o.sizeX) + "x" + std::to_string(o.sizeY) + ")");
    for (size_t i = 0; i < this->size(); i++)
        this->_data[i] += o._data[i];
    return *this;
}
Matrix operator-(const Matrix& a, const Matrix& b)
{
    Matrix temp = a;
    temp -= b;
    return temp;
}
Matrix& Matrix::operator-=(const Matrix& o)
{
    if (sizeX != o.sizeX || sizeY != o.sizeY)
        throw std::domain_error("Error on matrix substraction. Matrix A (" + std::to_string(sizeX) + "x" + std::to_string(sizeY) + ") does not match Matrix B (" + std::to_string(o.sizeX) + "x" + std::to_string(o.sizeY) + ")");
    for (size_t i = 0; i < this->size(); i++)
        this->_data[i] -= o._data[i];
    return *this;
}

Matrix Matrix::matprod(const Matrix &A, const Matrix &B)
{
    return A.product(B);
}

Matrix Matrix::toHomogeneous() const
{
    if (cols() != 3 || rows() != 3)
        throw std::invalid_argument("Homogeneous conversion requires 3x3 matrix");
    Matrix res(4, 4);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            res(i, j) = this->unsafe(i, j);
        }
    }
    res(3, 3) = 1;
    return res;
}

std::pair<Matrix, Matrix> Matrix::gramSchmidtQR(const Matrix &A)
{
    validate(A);
    const size_t m = A.rows(), n = A.cols();
    if (m < n)
        throw std::invalid_argument("QR requires rows >= columns");
    Matrix Q(n, m), R(n, n);
    std::vector<double> v(m);
    for (size_t k = 0; k < n; ++k) {
        double original = 0;
        for (size_t r = 0; r < m; ++r) {
            v[r] = A.unsafe(k, r);
            original += v[r] * v[r];
        }
        // Modified Gram-Schmidt with reorthogonalization for nearly parallel columns.
        for (int pass = 0; pass < 2; ++pass)
            for (size_t j = 0; j < k; ++j) {
                double dot = 0;
                for (size_t r = 0; r < m; ++r)
                    dot += Q.unsafe(j, r) * v[r];
                R.unsafe(k, j) += float(dot);
                for (size_t r = 0; r < m; ++r)
                    v[r] -= dot * Q.unsafe(j, r);
            }
        double norm = 0;
        for (double x : v)
            norm += x * x;
        norm = std::sqrt(norm);
        if (!std::isfinite(norm) || norm <= std::numeric_limits<float>::epsilon() * std::sqrt(original))
            throw std::domain_error("Rank deficient matrix");
        R.unsafe(k, k) = float(norm);
        for (size_t r = 0; r < m; ++r)
            Q.unsafe(k, r) = float(v[r] / norm);
    }
    return {Q, R};
}

std::vector<float> Matrix::backSubstitution(const Matrix &R, const std::vector<float> &b)
{
    square(R);
    size_t n = R.rows();
    if (b.size() != n)
        throw std::invalid_argument("Right-hand side size mismatch");
    std::vector<float> x(n);
    for (size_t i = n; i-- > 0;) {
        double value = b[i];
        for (size_t j = i + 1; j < n; ++j)
            value -= R.unsafe(j, i) * double(x[j]);
        double d = R.unsafe(i, i);
        if (d == 0 || !std::isfinite(d))
            throw std::domain_error("Singular triangular system");
        x[i] = float(value / d);
    }
    return x;
}

std::vector<float> Matrix::solve(const Matrix &A, const Matrix &b)
{
    validate(b);
    if (b.cols() != 1 || b.rows() != A.rows())
        throw std::invalid_argument("Right-hand side must be a column matching A rows");
    return solve(A, b._data);
}

std::vector<float> Matrix::solve(const Matrix &A, const std::vector<float> &b)
{
    if (b.size() != A.rows())
        throw std::invalid_argument("Right-hand side size mismatch");
    for (float value : b)
        if (!std::isfinite(value))
            throw std::invalid_argument("Nonfinite right-hand side");
    auto [Q, R] = gramSchmidtQR(A);
    std::vector<float> rhs(A.cols());
    for (size_t c = 0; c < A.cols(); ++c) {
        double sum = 0;
        for (size_t r = 0; r < A.rows(); ++r)
            sum += Q.unsafe(c, r) * double(b[r]);
        rhs[c] = float(sum);
    }
    return backSubstitution(R, rhs);
}

Matrix operator*(const Matrix& a, const Matrix& o)
{
    return a.product(o);
}
Matrix& Matrix::operator*=(const Matrix& o)
{
    *this = Matrix::matprod(*this, o);
    return *this;
}
Matrix Matrix::operator*(float o)
{
    Matrix m = *this;
    m *= o;
    return m;
}
Matrix& Matrix::operator*=(float o)
{
    for (size_t i = 0; i < this->size(); i++)
        this->_data[i] *= o;
    return *this;
}
Matrix Matrix::operator/(float o)
{
    Matrix m(*this);
    m /= o;
    return m;
}
Matrix& Matrix::operator/=(float o)
{
    if (o == 0 || !std::isfinite(o))
        throw std::domain_error("Invalid matrix divisor");
    for (size_t i = 0; i < this->size(); i++)
        this->_data[i] /= o;
    return *this;
}
Matrix Matrix::operator+(float o)
{
    Matrix m(*this);
    m += o;
    return m;
}
Matrix& Matrix::operator+=(float o)
{
    for (size_t i = 0; i < this->size(); i++)
        this->_data[i] += o;
    return *this;
}
Matrix Matrix::operator-(float o)
{
    Matrix m(*this);
    m -= o;
    return m;
}
Matrix& Matrix::operator-=(float o)
{
    for (size_t i = 0; i < this->size(); i++)
        this->_data[i] -= o;
    return *this;
}


std::ostream& operator<<(std::ostream& io, const Matrix& m) {
    io << "Matrix (" << m.cols() << "x" << m.rows() << ") :\n";
    for (size_t row = 0; row < m.rows(); row ++)
    {
        for (size_t col = 0; col < m.cols(); col++) {
            io << m(col, row) << "\t";
        }
        io << "\n";
    }
    return io;
}

std::ostream& operator<<(std::ostream& io, Matrix* m) {
    io << m->toString();
    return io;
}

Matrix Matrix::operator/(const Matrix &other) {
    return product(other.inverse());
}
Matrix &Matrix::operator/=(const Matrix &other) {
    return *this = *this / other;
}
