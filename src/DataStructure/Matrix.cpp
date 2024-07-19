#include "DataStructure/Matrix.h"


#include <cmath>
#include <sstream>

Matrix::Matrix()
    : Matrix(1, 1, 0.0)
{

}

Matrix::Matrix(int n, int m, float default_value)
{
    for(int i = 0; i < n; i++) {
        this->push_back(std::vector<float>(m, default_value));
    }
}
Matrix::Matrix(int n, int m, float* data){
    for(int i = 0; i < n; i++) {
        this->push_back(std::vector<float>());
        for (int j = 0; j < m; j++)
            (*this)[i].push_back(data[i + j*n]);
    }
}
Matrix::Matrix(int n, int m, float** data){

    for(int i = 0; i < n; i++) {
        this->push_back(std::vector<float>());
        for (int j = 0; j < m; j++)
            (*this)[i].push_back(data[i][j]);
    }
}
Matrix::Matrix(const std::vector<std::vector<float> > &data)
    : Matrix(data.size(), data[0].size())
{
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0].size(); j++) {
            (*this)[i][j] = data[i][j];
        }
    }
//    (*this) = data;
}

float Matrix::det() const
{
    float determinant = 0;
    if (this->size() == 1) {
       return (*this)[0][0];
    }
    if (this->size() == 2) {
       return ((*this)[0][0] * (*this)[1][1]) - ((*this)[0][1] * (*this)[1][0]);
    }
    float sign = 1.0;
    for (size_t i = 0; i < this->size(); i++) {
       Matrix temp = this->submatrix(0, i);
       determinant += sign * (*this)[0][i] * temp.det();
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
    for (size_t row = 0; row < (*this)[0].size(); row++)
    {
        for (size_t col = 0; col < this->size(); col++)
        {
            temp[row][col] = this->submatrix(row, col).det();
        }
    }
    return temp;
}
Matrix Matrix::submatrix(size_t rowToIgnore, size_t colToIgnore) const
{
    Matrix temp(this->size() - 1, this->size() - 1);
    int i = 0, j = 0;
    // filling the sub matrix
    for (size_t row = 0; row < this->size(); row++) {
       for (size_t col = 0; col < this->size(); col++) {
          // skipping if the current row or column is not equal to the current
          // element row and column
          if (row != rowToIgnore && col != colToIgnore) {
             temp[i][j++] = (*this)[row][col];
             if (j == int(this->size() - 1)) {
                j = 0;
                i++;
             }
          }
       }
    }
    return temp;
}

Matrix Matrix::transpose() const
{
    Matrix temp((*this)[0].size(), this->size());
    for(size_t row = 0; row < (*this)[0].size(); row++)
        for(size_t col = 0; col < this->size(); col++)
            temp[row][col] = (*this)[col][row];
    return temp;
}

Matrix Matrix::product(const Matrix& mat) const
{
    int m = this->cols();
    int l = this->rows();
    int m2 = mat.rows();
    int n = mat.cols();
    if (m != m2) {
        std::ostringstream oss;
        oss << "Error on product between matrices. A defined as " << m << "x" << l << " and B defined as " << n << "x" << m2 << ":\n";
        oss << this->toString() << "\n*\n" << mat.toString();
        throw std::domain_error(oss.str());
    }
    Matrix temp(l, n);
    for (size_t row = 0; row < l; row++) {
        for (size_t col = 0; col < n; col++)
        {
            for (size_t k = 0; k < m; k++) {
                temp[row][col] += (*this)[row][k] * mat[k][col];
            }
        }
    }
    return temp;
}

float Matrix::trace() const
{
    float trace = 0.0;
    for (size_t i = 0; i < this->size(); i++)
        trace += (*this)[i][i];
    return trace;
}

Matrix Matrix::col(int colIndex) const
{
    Matrix res(1, this->rows());
    for (int i = 0; i < this->rows(); i++) {
        res[0][i] = (*this)[colIndex][i];
    }
    return res;
}

Matrix Matrix::row(int rowIndex) const
{
    Matrix res(this->cols(), 1);
    for (int i = 0; i < this->cols(); i++) {
        res[i][0] = (*this)[i][rowIndex];
    }
    return res;
}

Matrix Matrix::abs() const
{
    Matrix res = *this;
    for (int i = 0; i < res.size(); i++) {
        for (int j = 0; j < res[i].size(); j++) {
            res[i][j] = std::abs(res[i][j]);
        }
    }
    return res;
}

float Matrix::maxCoeff() const
{
    float max = std::numeric_limits<float>::lowest();

    for (int i = 0; i < (*this).size(); i++) {
        for (int j = 0; j < (*this)[i].size(); j++) {
            max = std::max(max, (*this)[i][j]);
        }
    }
    return max;
}

Matrix Matrix::leftCols(int nbCols) const
{
    Matrix res(this->rows(), nbCols);
    for (int i = 0; i < res.size(); i++) {
        for (int j = 0; j < res[i].size(); j++) {
            res[i][j] = (*this)[i][j];
        }
    }
    return res;
}

Matrix Matrix::rightCols(int nbCols) const
{
    Matrix res(this->cols() - nbCols, this->rows());
    for (int i = 0; i < res.size(); i++) {
        for (int j = 0; j < res[i].size(); j++) {
            res[i][j] = (*this)[i + nbCols][j];
        }
    }
    return res;
}

Matrix Matrix::identity(int size)
{
    Matrix m(size, size);
    for (int i = 0; i < size; i++) {
        m[i][i] = 1;
    }
    return m;
}

std::string Matrix::displayValues() const
{
    std::string txt = "";
    for (size_t col = 0; col < this->size(); col ++) {
        for (size_t row = 0; row < this[0].size(); row ++) {
            txt += std::to_string(int((*this)[row][col] * 1000)/1000.) + "\t";
        }
        txt += "\n";
    }
    return txt;
}

std::string Matrix::toString() const
{
    return "Matrix (" + std::to_string(this->size()) + "x" + std::to_string(this[0].size()) + ") :\n" + this->displayValues();
    //    return txt;
}

std::string Matrix::displayValuesOneLine() const
{
    std::string txt = "";
    for (size_t col = 0; col < this->size(); col ++) {
        for (size_t row = 0; row < this[0].size(); row ++) {
            txt += std::to_string(int((*this)[row][col] * 1000)/1000.) + " ";
        }
    }
    return txt;
}

float* Matrix::toArray(bool rowByRow) const
{
    const int nbCols = cols();
    const int nbRows = rows();

    float* values = new float(nbRows * nbCols);
    for (size_t r = 0; r < nbRows; r++) {
        for (size_t c = 0; c < nbCols; c++) {
            if (rowByRow)
                values[r * nbCols + c] = (*this)[r][c]; // Row by row
            else
                values[c * nbRows + r] = (*this)[r][c]; // Column by column
        }
    }
    return values;
}

std::vector<float> Matrix::toStdVector(bool rowByRow) const
{
    const int nbCols = cols();
    const int nbRows = rows();

    std::vector<float> values(nbRows * nbCols);
    for (size_t r = 0; r < nbRows; r++) {
        for (size_t c = 0; c < nbCols; c++) {
            if (rowByRow)
                values[r * nbCols + c] = (*this)[r][c]; // Row by row
            else
                values[c * nbRows + r] = (*this)[r][c]; // Column by column
        }
    }
    return values;
}


Matrix operator+(const Matrix& a, const Matrix& b)
{
    Matrix temp = a;
    temp += b;
    return temp;
}
Matrix& Matrix::operator+=(const Matrix& o)
{
    for(size_t row = 0; row < o.size(); row++)
        for(size_t col = 0; col < o[0].size(); col++)
            (*this)[row][col] += o[row][col];
    return (*this);
}
Matrix operator-(const Matrix& a, const Matrix& b)
{
    Matrix temp = a;
    temp -= b;
    return temp;
}
Matrix& Matrix::operator-=(const Matrix& o)
{
    for(size_t row = 0; row < o.size(); row++)
        for(size_t col = 0; col < o[0].size(); col++)
            (*this)[row][col] -= o[row][col];
    return (*this);
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
            res[i][j] = (*this)[i][j];
        }
    }
    res[3][3] = 1;
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
            Q[i][k] = A[i][k];
        }

        for (size_t j = 0; j < k; ++j) {
            float dot_product = 0;
            for (size_t i = 0; i < rows; ++i) {
                dot_product += Q[i][j] * A[i][k];
            }
            R[j][k] = dot_product;

            for (size_t i = 0; i < rows; ++i) {
                Q[i][k] -= dot_product * Q[i][j];
            }
        }

        float norm = 0;
        for (size_t i = 0; i < rows; ++i) {
            norm += Q[i][k] * Q[i][k];
        }
        norm = std::sqrt(norm);
        R[k][k] = norm;

        for (size_t i = 0; i < rows; ++i) {
            Q[i][k] /= norm;
        }
    }

    return {Q, R};
}

std::vector<float> Matrix::backSubstitution(const Matrix &R, const vector<float> &b)
{
    size_t n = R.rows();
    vector<float> x(n);

    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (size_t j = i + 1; j < n; ++j) {
            x[i] -= R[i][j] * x[j];
        }
        x[i] /= R[i][i];
    }

    return x;
}

std::vector<float> Matrix::solve(const Matrix &A, const Matrix &b)
{
    auto [Q, R] = Matrix::gramSchmidtQR(A);
    Matrix Qt = Q.transpose();
    vector<float> Qt_b(b.size());

    for (size_t i = 0; i < Qt.rows(); ++i) {
        for (size_t j = 0; j < Qt.cols(); ++j) {
            Qt_b[i] += Qt[i][j] * b[j][0];
        }
    }

    return backSubstitution(R, Qt_b);
}

std::vector<float> Matrix::solve(const Matrix &A, const std::vector<float> &b)
{
    Matrix B(b.size(), 1);
    for (int i = 0; i < b.size(); i++) {
        B[i][0] = b[i];
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
    *this = Matrix::matprod(*this, o);/*
    for(size_t row = 0; row < o.size(); row++)
        for(size_t col = 0; col < o[0].size(); col++)
            (*this)[row][col] *= o[row][col];*/
    return (*this);
}
Matrix Matrix::operator*(float o)
{
    Matrix m(*this);
    m *= o;
    return m;
}
Matrix& Matrix::operator*=(float o)
{
    for(size_t row = 0; row < this->size(); row++)
        for(size_t col = 0; col < (*this)[0].size(); col++)
            (*this)[row][col] *= o;
    return (*this);
}
Matrix Matrix::operator/(float o)
{
    Matrix m(*this);
    m /= o;
    return m;
}
Matrix& Matrix::operator/=(float o)
{
    for(size_t row = 0; row < this->size(); row++)
        for(size_t col = 0; col < (*this)[0].size(); col++)
            (*this)[row][col] /= o;
    return (*this);
}
Matrix Matrix::operator+(float o)
{
    Matrix m(*this);
    m += o;
    return m;
}
Matrix& Matrix::operator+=(float o)
{
    for(size_t row = 0; row < this->size(); row++)
        for(size_t col = 0; col < (*this)[0].size(); col++)
            (*this)[row][col] += o;
    return (*this);
}
Matrix Matrix::operator-(float o)
{
    Matrix m(*this);
    m -= o;
    return m;
}
Matrix& Matrix::operator-=(float o)
{
    for(size_t row = 0; row < this->size(); row++)
        for(size_t col = 0; col < (*this)[0].size(); col++)
            (*this)[row][col] -= o;
    return (*this);
}


std::ostream& operator<<(std::ostream& io, const Matrix& m) {
    io << "Matrix (" << m.size() << "x" << m[0].size() << ") :\n";
    for (size_t row = 0; row < m.size(); row ++)
    {
        for (size_t col = 0; col < m[0].size(); col ++) {
            io << int(m[row][col] * 1000)/1000. << "\t";
        }
        io << "\n";
    }
    return io;
}

std::ostream& operator<<(std::ostream& io, Matrix* m) {
    io << m->toString();
    return io;
}
