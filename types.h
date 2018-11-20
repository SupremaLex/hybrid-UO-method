#ifndef TYPES
#define TYPES
#include <cmath>
#include <vector>
#include <deque>
#include <string>
//#define ARMA_DONT_USE_OPENMP
#include <armadillo>

using std::vector;
using std::deque;
using std::string;

typedef arma::Mat<double> Matrix;
typedef arma::Col<double> Vector;
typedef double (*Function)(const Vector&);
#endif // TYPES

