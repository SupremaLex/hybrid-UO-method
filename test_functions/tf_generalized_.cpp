#include "test_function_collection.h"
namespace gen = test_functions::generalized;


double gen::gen_rosenbrock(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n - 1, f)
      100 * pow(x_i(i+1) - pow(x_i(i), 2), 2) + pow(1 - x_i(i), 2);
  return f;
}


double gen::gen_white_holst(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n - 1, f)
      100 * pow(x_i(i+1) - pow(x_i(i), 3), 2) + pow(1 - x_i(i), 2);
  return f;
}


double gen::gen_PSC1(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n - 1, f)
      pow(pow(x_i(i), 2) + pow(x_i(i+1), 2) + x_i(i) * x_i(i+1), 2) +
      pow(sin(x_i(i)), 2) + pow(cos(x_i(i)), 2);
  return f;
}
