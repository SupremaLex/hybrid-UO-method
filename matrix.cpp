#include "matrix.h"
#include <exception>

Matrix::Matrix(int n, int m) {
  _matrix_repr = vector< vector<T>(m) > (n);
  _rows = n;
  _colomns = m;
}


Matrix::Matrix(vector<vector<T> > &m) {
   _matrix_repr = m;
   _rows = m.size();
   _colomns = m[0].size();
}


Matrix& Matrix::operator +(const Matrix& m) {
  if (m.size() != _rows && m.at(0).size() != _colomns)
    throw std::exception("Matrix dimensions not equal");
  for (int i = 0; i < _rows; i++)
    for (int j = 0; j < _colomns; j++)
      _matrix_repr[i][j] += m[i][j];
  return *this;
}


Matrix& Matrix::operator -(const Matrix& m) {
  if (m.size() != _rows && m.at(0).size() != _colomns)
    throw std::exception("Matrix dimensions not equal");
  for (int i = 0; i < _rows; i++)
    for (int j = 0; j < _colomns; j++)
      _matrix_repr[i][j] -= m[i][j];
  return *this;
}


Matrix& Matrix::operator *(const Matrix& m) {
  if (m.size() != _colomns)
    throw std::exception("Matrix dimensions wrong");
  T sum;
  for (int i = 0; i < _rows; i++)
    for (int j = 0; j < _colomns; j++) {
      sum = 0;
      for (int k = 0; k < _colomns; k++)
        sum += _matrix_repr[i][k] * m[k][j];
      _matrix_repr[i][j] = sum;
    }

   return *this;
}

ostream& Matrix::operator <<(ostream& os) {
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _colomns; j++)
      os << _matrix_repr[i][j];
    os << std::endl;
  }
  return os;
}
