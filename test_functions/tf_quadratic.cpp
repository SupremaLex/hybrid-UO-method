#include "test_function_collection.h"
namespace quad = test_functions::quadratic;


double quad::quadratic_1(const Vector& x) {
  double f = 0, f1 = 0;
  int n = x.size();
  sum(1, n, f)
      x_i(i);
  sum(1, n, f1)
      i * pow(x_i(i), 2) / 100;
  return f * f + f1;
}


double quad::quadratic_2(const Vector& x) {
  int n = x.size();
  double f = 0;
  sum(1, n, f)
      i * pow(pow(x_i(i) - 1, 2) - 1, 2);
  return 0.5 * f + x_i(n);
}


double quad::perturbed_quadratic(const Vector& x) {
  double f = 0, f1 = 0;
  int n = x.size();
  sum(1, n, f)
      i * pow(x_i(i), 2);
  sum(1, n, f1)
      x_i(i);
  return f + f1 * f1 / 100;
}


double quad::perturbed_quadratic_diagonal(const Vector& x) {
  double f = 0, f1 = 0;
  int n = x.size();
  sum(1, n, f)
      pow(x_i(i), 2) * i;
  sum(1, n, f1)
      x_i(i);
  return f / 100 + f1 * f1;
}


double quad::ext_quadratic_penalty_1(const Vector& x) {
  double f = 0, f1 = 0;
  int n = x.size();
  sum(1, n - 1, f)
      pow(pow(x_i(i), 2) - 2, 2);
  sum(1, n, f1)
      pow(x_i(i), 2);
  return f + (f1 - 0.5) * (f1 - 0.5);
}


double quad::ext_quadratic_penalty_2(const Vector& x) {
  double f = 0, f1 = 0;
  int n = x.size();
  sum(1, n - 1, f)
      pow(pow(x_i(i), 2) - sin(x_i(i)), 2);
  sum(1, n, f1)
      pow(x_i(i), 2);
  return f + (f1 - 100) * (f1 - 100);
}


double quad::ext_quadratic_exponential_1(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n / 2, f)
      pow(exp(x_i(2*i-1) - x_i(2*i) - 5), 2) +
      pow(x_i(2*i-1) - x_i(2*i), 2) * pow(x_i(2*i-1) - x_i(2*i) - 11, 2);
  return f;
}


double quad::partial_perturbed_quadratic(const Vector& x) {
  double f = pow(x_i(1), 2);
  int n = x.size();
  double tmp = 0;
  for (int i = 1; i < n + 1; i++) {
    tmp += x_i(i);
    f += i * pow(x_i(i), 2) + tmp / 100;
  }
  return f;
}


double quad::almost_perturbed_quadratic(const Vector& x) {
  int n = x.size();
  double f = pow(x_i(1) + x_i(n), 2) / 100;
  sum(1, n, f)
      i * pow(x_i(i), 2);
  return f;
}


double quad::perturbed_tridiagonal_quadratic(const Vector& x) {
  double f = pow(x_i(1), 2);
  int n = x.size();
  sum(2, n - 1, f)
      i * pow(x_i(i), 2) + pow(x_i(i-1) + x_i(i) + x_i(i+1), 2);
  return f;
}


double quad::gen_quadratic(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n - 1, f)
      pow(x_i(i), 2) + pow(x_i(i+1) + pow(x_i(i), 2), 2);
  return f;
}
