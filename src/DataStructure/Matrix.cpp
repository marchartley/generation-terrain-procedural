#include "DataStructure/Matrix.h"

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
Matrix::Matrix(std::vector<std::vector<float>> data)
{
    (*this) = data;
}

float Matrix::det()
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
Matrix Matrix::submatrix(size_t rowToIgnore, size_t colToIgnore)
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

Matrix Matrix::transpose()
{
    Matrix temp(this[0].size(), this->size());
    for(size_t row = 0; row < this->size(); row++)
        for(size_t col = 0; col < this->size(); col++)
            temp[row][col] = (*this)[col][row];
    return temp;
}

Matrix Matrix::product(Matrix m)
{
    Matrix temp((*this)[0].size(), m[0].size());
    for (size_t row = 0; row < (*this)[0].size(); row++) {
        for (size_t col = 0; col < m[0].size(); col++)
        {
            for (size_t k = 0; k < m.size(); k++) {
                temp[row][col] += (*this)[row][k] * m[k][col];
            }
        }
    }
    return temp;
}

float Matrix::trace()
{
    float trace = 0.0;
    for (size_t i = 0; i < this->size(); i++)
        trace += (*this)[i][i];
    return trace;
}

std::string Matrix::displayValues()
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

std::string Matrix::toString()
{
    return "Matrix (" + std::to_string(this->size()) + "x" + std::to_string(this[0].size()) + ") :\n" + this->displayValues();
    //    return txt;
}

std::string Matrix::displayValuesOneLine()
{
    std::string txt = "";
    for (size_t col = 0; col < this->size(); col ++) {
        for (size_t row = 0; row < this[0].size(); row ++) {
            txt += std::to_string(int((*this)[row][col] * 1000)/1000.) + " ";
        }
    }
    return txt;
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
