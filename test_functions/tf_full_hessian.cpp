#include "test_function_collection.h"
namespace fh = test_functions::full_hessian;

double fh::full_hessian_1(const Vector& x) {
  double f = pow(x_i(1) - 3, 2);
  int n = x.size();
  double tmp = x_i(1);
  for (int i = 2; i < n; i++) {
      tmp += x_i(i);
      f += pow(x_i(1) - 3 + 2 * pow(tmp, 2), 2);
  }
  return f;
}

double fh::full_hessian_2(const Vector& x) {
  double f = pow(x_i(1) - 5, 2);
  int n = x.size();
  double tmp = x_i(1);
  for (int i = 2; i < n; i++) {
      tmp += x_i(i);
      f += pow(tmp - 1, 2);
  }
  return f;
}

double fh::full_hessian_3(const Vector& x) {
  double f = 0, f1 = 0;
  int n = x.size();
  sum(1, n, f1)
      x_i(i);
  sum(1, n, f)
      x_i(i) * exp(x_i(i)) - 2 * x_i(i) - pow(x_i(i), 2);
  return f + f1 * f1;
}
