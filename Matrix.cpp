#include "Matrix.h"

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
Matrix::Matrix(int n, int m, float** data){

    for(int i = 0; i < n; i++) {
        this->push_back(std::vector<float>());
        for (int j = 0; j < m; j++)
            (*this)[i].push_back(data[i][j]);
    }
}
Matrix::Matrix(std::vector<std::vector<float>> data)
{
    (*this) = data;
}

float Matrix::det()
{
    int determinant = 0;
    if (this->size() == 1) {
       return (*this)[0][0];
    }
    if (this->size() == 2) {
       return ((*this)[0][0] * (*this)[1][1]) - ((*this)[0][1] * (*this)[1][0]);
    }
    int sign = 1;
    for (int i = 0; i < this->size(); i++) {
       Matrix temp = this->submatrix(0, i);
       determinant += sign * (*this)[0][i] * temp.det();
       sign = -sign;
    }
    return determinant;
}
Matrix Matrix::adj()
{
    return this->cofactors().transpose();
}
Matrix Matrix::inverse()
{
    return this->adj() / this->det();
}
Matrix Matrix::cofactors()
{
    int a = 0;
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
Matrix Matrix::submatrix(int rowToIgnore, int colToIgnore)
{
    Matrix temp(this->size() - 1, this->size() - 1);
    int i = 0, j = 0;
    // filling the sub matrix
    for (int row = 0; row < this->size(); row++) {
       for (int col = 0; col < this->size(); col++) {
          // skipping if the current row or column is not equal to the current
          // element row and column
          if (row != rowToIgnore && col != colToIgnore) {
             temp[i][j++] = (*this)[row][col];
             if (j == this->size() - 1) {
                j = 0;
                i++;
             }
          }
       }
    }
    return temp;
}

Matrix Matrix::transpose()
{
    Matrix temp(this[0].size(), this->size());
    for(size_t row = 0; row < this->size(); row++)
        for(size_t col = 0; col < this->size(); col++)
            temp[row][col] = (*this)[col][row];
    return temp;
}


Matrix operator+(Matrix a, Matrix b)
{
    a += b;
    return a;
}
Matrix& Matrix::operator+=(const Matrix& o)
{
    for(size_t row = 0; row < o.size(); row++)
        for(size_t col = 0; col < o[0].size(); col++)
            (*this)[row][col] += o[row][col];
    return (*this);
}
Matrix operator-(Matrix a, const Matrix& b)
{
    a -= b;
    return b;
}
Matrix& Matrix::operator-=(const Matrix& o)
{
    for(size_t row = 0; row < o.size(); row++)
        for(size_t col = 0; col < o[0].size(); col++)
            (*this)[row][col] -= o[row][col];
    return (*this);
}
Matrix operator*(Matrix a, const Matrix& o)
{
    a *= o;
    return a;
}
Matrix& Matrix::operator*=(const Matrix& o)
{
    for(size_t row = 0; row < o.size(); row++)
        for(size_t col = 0; col < o[0].size(); col++)
            (*this)[row][col] *= o[row][col];
    return (*this);
}
Matrix Matrix::operator/(const Matrix& o)
{
    return (*this) / o;
}
Matrix& Matrix::operator/=(const Matrix& o)
{
    for(size_t row = 0; row < o.size(); row++)
        for(size_t col = 0; col < o[0].size(); col++)
            (*this)[row][col] /= o[row][col];
    return (*this);
}
Matrix Matrix::operator*(float o)
{
    return (*this) * o;
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
    return (*this) / o;
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
    return (*this) + o;
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
    return (*this) - o;
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
    for (size_t row = 0; row < m[0].size(); row ++)
    {
        for (size_t col = 0; col < m.size(); col ++) {
            io << m[row][col] << "\t";
        }
        io << "\n";
    }
    return io;
}

std::ostream& operator<<(std::ostream& io, Matrix* m) {
    io << "Matrix (" << m->size() << "x" << m[0].size() << ") :\n";
    for (size_t row = 0; row < m[0].size(); row ++)
    {
        for (size_t col = 0; col < m->size(); col ++) {
            io << (*m)[row][col] << "\t";
        }
        io << "\n";
    }
    return io;
}
