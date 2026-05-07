#include "DataStructure/Matrix.h"


#include <cmath>
#include <sstream>

Matrix::Matrix()
    : Matrix(1, 1, 0.0)
{}

Matrix::Matrix(size_t n, size_t m, float default_value)
{
    this->sizeX = n;
    this->sizeY = m;
    this->_data = std::vector<float>(sizeX * sizeY, default_value);
}
Matrix::Matrix(size_t n, size_t m, const float* data){
    this->sizeX = n;
    this->sizeY = m;
    this->_data = std::vector<float>(sizeX * sizeY);
    this->_data.assign(data, data + sizeX * sizeY);
}

Matrix::Matrix(size_t n, size_t m, const std::vector<float> &data)
    : Matrix(n, m)
{
    this->_data = data;
}
Matrix::Matrix(const std::vector<std::vector<float> > &data)
    : Matrix(data[0].size(), data.size())
{
    for (size_t i = 0; i < data.size(); i++) {
        for (size_t j = 0; j < data[0].size(); j++) {
            this->unsafe(j, i) = data[i][j];
        }
    }
}

float Matrix::det() const
{
    float determinant = 0;
    if (cols() == 1) {
        return this->unsafe(0, 0);
    }
    if (cols() == 2) {
        return (this->unsafe(0, 0) * this->unsafe(1, 1)) - (this->unsafe(0, 1) * this->unsafe(1, 0));
    }
    float sign = 1.0;
    for (size_t i = 0; i < cols(); i++) {
        Matrix temp = this->submatrix(0, i);
        determinant += sign * this->unsafe(0, i) * temp.det();
        sign = -sign;
    }
    return determinant;
}
Matrix Matrix::adj() const
{
    return this->cofactors().transpose();
}
Matrix Matrix::inverse() const
{
    return this->adj() / this->det();
}
Matrix Matrix::cofactors() const
{
    Matrix temp(*this);
    // Looping for each element of the matrix
    for (size_t row = 0; row < rows(); row++)
    {
        for (size_t col = 0; col < cols(); col++)
        {
            temp.unsafe(col, row) = this->submatrix(col, row).det();
        }
    }
    return temp;
}
Matrix Matrix::submatrix(size_t rowToIgnore, size_t colToIgnore) const
{
    Matrix temp(cols() - 1, rows() - 1);
    size_t i = 0, j = 0;
    // filling the sub matrix
    for (size_t row = 0; row < rows(); row++) {
        for (size_t col = 0; col < cols(); col++) {
            // skipping if the current row or column is not equal to the current
            // element row and column
            if (row == rowToIgnore || col == colToIgnore) continue;

            temp.unsafe(i, j++) = this->unsafe(col, row);
            if ((int)j == int(cols() - 1)) {
                j = 0;
                i++;
            }
        }
    }
    return temp;
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
                temp(col, row) += a * mat.unsafe(col, k);
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
    Matrix res(1, this->rows());
    for (size_t i = 0; i < this->rows(); i++) {
        res._data[i] = this->unsafe(colIndex, i);
    }
    return res;
}

Matrix Matrix::row(size_t rowIndex) const
{
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
    Matrix res(this->rows(), nbCols);
    for (size_t i = 0; i < res.cols(); i++) {
        for (size_t j = 0; j < res.rows(); j++) {
            res(i, j) = this->unsafe(i, j);
        }
    }
    return res;
}

Matrix Matrix::rightCols(size_t nbCols) const
{
    Matrix res(this->cols() - nbCols, this->rows());
    for (size_t i = 0; i < res.cols(); i++) {
        for (size_t j = 0; j < res.rows(); j++) {
            res(i, j) = this->unsafe(i + nbCols, j);
        }
    }
    return res;
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
            txt += std::to_string(int(this->unsafe(col, row) * 1000)/1000.) + "\t";
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
            txt += std::to_string(int(this->unsafe(col, row) * 1000)/1000.) + " ";
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
    size_t rows = A.rows();
    size_t cols = A.cols();
    Matrix Q(rows, cols);
    Matrix R(cols, cols);

    for (size_t k = 0; k < cols; ++k) {
        for (size_t i = 0; i < rows; ++i) {
            Q(i, k) = A(i, k);
        }

        for (size_t j = 0; j < k; ++j) {
            float dot_product = 0;
            for (size_t i = 0; i < rows; ++i) {
                dot_product += Q(i, j) * A(i, k);
            }
            R(j, k) = dot_product;

            for (size_t i = 0; i < rows; ++i) {
                Q(i, k) -= dot_product * Q(i, j);
            }
        }

        float norm = 0;
        for (size_t i = 0; i < rows; ++i) {
            norm += Q(i, k) * Q(i, k);
        }
        norm = std::sqrt(norm);
        R(k, k) = norm;

        for (size_t i = 0; i < rows; ++i) {
            Q(i, k) /= norm;
        }
    }

    return {Q, R};
}

std::vector<float> Matrix::backSubstitution(const Matrix &R, const std::vector<float> &b)
{
    size_t n = R.rows();
    std::vector<float> x(n);

    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (size_t j = i + 1; j < n; ++j) {
            x[i] -= R(i, j) * x[j];
        }
        x[i] /= R(i, i);
    }

    return x;
}

std::vector<float> Matrix::solve(const Matrix &A, const Matrix &b)
{
    auto [Q, R] = Matrix::gramSchmidtQR(A);
    Matrix Qt = Q.transpose();
    std::vector<float> Qt_b(b.rows());

    for (size_t i = 0; i < Qt.rows(); ++i) {
        for (size_t j = 0; j < Qt.cols(); ++j) {
            Qt_b[i] += Qt(i, j) * b(j, 0);
        }
    }

    return backSubstitution(R, Qt_b);
}

std::vector<float> Matrix::solve(const Matrix &A, const std::vector<float> &b)
{
    Matrix B(b.size(), 1);
    for (size_t i = 0; i < b.size(); i++) {
        B(i, 0) = b[i];
    }
    return Matrix::solve(A, B);
}

Matrix operator*(const Matrix& a, const Matrix& o)
{
    Matrix temp = a;
    temp *= o;
    return temp;
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
        for (size_t col = 0; col < m.rows(); col ++) {
            io << int(m(col, row) * 1000)/1000. << "\t";
        }
        io << "\n";
    }
    return io;
}

std::ostream& operator<<(std::ostream& io, Matrix* m) {
    io << m->toString();
    return io;
}
