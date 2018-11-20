#include "test_function_collection.h"
namespace diag = test_functions::diagonal;


double diag::diagonal_1(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n, f)
      exp(x_i(i)) - i * x_i(i);
  return f;
}


double diag::diagonal_2(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n, f)
      exp(x_i(i)) - x_i(i) / i;
  return f;
}


double diag::diagonal_3(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n, f)
      exp(x_i(i)) - i * sin(x_i(i));
  return f;
}


double diag::diagonal_4(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n / 2, f)
      0.5 * (pow(x_i(2*i-1), 2) + 100 * pow(x_i(2*i), 2));
  return f;
}


double diag::diagonal_5(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n, f)
      log(exp(x_i(i)) + exp(-x_i(i)));
  return f;
}


double diag::diagonal_6(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n, f)
      exp(x_i(i)) - 1 + x_i(i);
  return f;
}


double diag::diagonal_7(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n, f)
      exp(x_i(i)) - 2 * x_i(i) - pow(x_i(i), 2);
  return f;
}


double diag::diagonal_8(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n, f)
      x_i(i) * exp(x_i(i)) - 2 * x_i(i) - pow(x_i(i), 2);
  return f;
}


double diag::diagonal_9(const Vector& x) {
  int n = x.size();
  double f = 10000 * x_i(n);
  sum(1, n, f)
      exp(x_i(i)) - i * x_i(i);
  return f;
}


double diag::ext_tridiagonal(const Vector& x) {
  int n = x.size();
  double f = 0;
  sum(1, n / 2, f)
      pow(x_i(2*i-1) + x_i(2*i) - 3, 2) +
      pow(x_i(2*i-1) - x_i(2*i) + 1, 4);
  return f;
}


double diag::ext_block_diagonal_1(const Vector& x) {
  int n = x.size();
  double f = 0;
  sum(1, n / 2, f)
      pow(pow(x_i(2*i-1), 2) + pow(x_i(2*i), 2) - 2, 2) +
      pow(exp(x_i(2*i-1) - 1) - x_i(2*i), 2);
  return f;
}


double diag::ext_tridiagonal_2(const Vector& x) {
  int n = x.size();
  double f = 0;
  sum(1, n - 1, f)
      pow(x_i(i) * x_i(i+1) - 1, 2) + 0.1 * (x_i(i) + 1) * (x_i(i+1) + 1);
  return f;
}


double diag::broyden_tridiagonal(const Vector& x) {
  int n = x.size();
  double f = pow(3 * x_i(1) - 2 *pow(x_i(1), 2), 2) + pow(3 * x_i(n) - 2 * pow(x_i(n), 2) - x_i(n-1) + 1, 2);
  sum(2, n - 1, f)
      pow(3 * x_i(i) - 2 * pow(x_i(i), 2) - x_i(i) - 2 * x_i(i+1) + 1, 2);
  return f;
}


double diag::gen_tridiagonal_1(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n - 1, f)
      pow(x_i(i) + x_i(i+1) - 3, 2) + pow(x_i(i) - x_i(i+1) + 1, 4);
  return f;
}


double diag::gen_tridiagonal_2(const Vector& x) {
  int n = x.size();
    double f = pow((5 - 3 *x_i(1) - pow(x_i(1), 2)) * x_i(1) - 3 * x_i(2) + 1, 2) +
        pow((5 - 3 *x_i(n) - pow(x_i(n), 2)) * x_i(n) - x_i(n-1) + 1, 2);
    sum(2, n - 1, f)
        pow((5 - 3 *x_i(i) - pow(x_i(i), 2)) * x_i(i) - x_i(i-1) - 3 * x_i(i+1) + 1, 2);
    return f;
}
