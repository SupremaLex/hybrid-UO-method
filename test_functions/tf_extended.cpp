#include "test_function_collection.h"
namespace ext = test_functions::extended;

double ext::ext_freudenstein_roth(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n / 2, f)
      pow(-13 + x_i(2*i-1) + ((5 - x_i(2*i))*x_i(2*i) - 2)*x_i(2*i), 2) +
      pow(-29 + x_i(2*i-1) + ((x_i(2*i) + 1) * x_i(2*i) - 14) * x_i(2*i), 2);
  return f;
}


double ext::ext_trigonometric(const Vector& x) {
  double f = 0;
  double tmp = 0;
  int n = x.size();
  sum(1, n, tmp)
      cos(x_i(i));
  sum(1, n, f)
      pow(n - tmp + i * (1 - cos(x_i(i))) - sin(x_i(i)), 2);
  return f;
}


double ext::ext_rosenbrock(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n / 2, f)
      100 * pow(x_i(2*i) - pow(x_i(2*i-1), 2), 2) + pow(1 - x_i(2*i-1), 2);
  return f;
}


double ext::ext_white_holst(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n / 2, f)
      100 * pow(x_i(2*i) - pow(x_i(2*i-1), 3), 2) + pow(1 - x_i(2*i-1), 2);
  return f;
}


double ext::ext_beale(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n / 2, f)
      pow(1.5 - x_i(2*i-1) * (1- x_i(2*i)), 2) +
      pow(2.25 - x_i(2*i-1) * (1 - pow(x_i(2*i), 2)), 2) +
      pow(2.625 - x_i(2*i-1) * (1 - pow(x_i(2*i), 3)), 2);
  return f;
}


double ext::ext_penalty(const Vector& x) {
  double f = 0;
  int n = x.size();
  double tmp = 0;
  sum(1, n, tmp)
      pow(x_i(i), 2);
  sum(1, n - 1, f)
      pow(x_i(i) - 1, 2);
  return f + (tmp - 0.25) * (tmp - 0.25);
}


double ext::ext_tet(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n / 2, f)
      exp(x_i(2*i-1) + 3 * x_i(2*i) - 0.1) +
      exp(x_i(2*i-1) - 3 * x_i(2*i) - 0.1) +
      exp(-x_i(2*i-1) - 0.1);
  return f;
}


double ext::ext_himmelblau(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n /2, f)
      pow(pow(x_i(2*i-1), 2) + x_i(2*i) - 11, 2) +
      pow(pow(x_i(2*i), 2) + x_i(2*i-1) - 7, 2);
  return f;
}


double ext::ext_PSC1(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n / 2, f)
      pow(pow(x_i(2*i-1), 2) + pow(x_i(2*i), 2) + x_i(2*i-1) * x_i(2*i), 2) +
      pow(sin(x_i(2*i-1)), 2) + pow(cos(x_i(2*i)), 2);
  return f;
}


double ext::ext_powell(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n / 4, f)
      pow(x_i(4*i-3) + 10 * x_i(4*i-2), 2) + 5 * pow(x_i(4*i-1) - x_i(4*i), 2) +
      pow(x_i(4*i-2) - 2 * x_i(4*i-1), 4) + 10 * pow(x_i(4*i-3) - x_i(4*i), 4);
  return f;
}


double ext::ext_maratos(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n / 2, f)
      x_i(2*i-1) + 100 * pow(pow(x_i(2*i-1), 2) + pow(x_i(2*i), 2) - 1, 2);
  return f;
}


double ext::ext_cliff(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n / 2, f)
      pow((x_i(2*i-1) - 3) / 100, 2) - (x_i(2*i-1) - x_i(2*i)) +
      exp(20 * (x_i(2*i-1) - x_i(2*i)));
  return f;
}


double ext::ext_wood(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n / 4, f)
      100 * pow(pow(x_i(4*i-3), 2) - x_i(4*i-2), 2) + pow(x_i(4*i-3) - 1, 2) +
      90 * pow(pow(x_i(4*i-1), 2) - x_i(4*i), 2) + pow(1 - x_i(4*i-1), 2) +
      10.1 * (pow(x_i(4*i-2) - 1, 2) + pow(x_i(4*i) - 1, 2)) +
      19.8 * (x_i(4*i-2) - 1) * (x_i(4*i) - 1);
  return f;
}


double ext::ext_hiebert(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n / 2, f)
      pow(x_i(2*i-1) - 10, 2) + pow(x_i(2*i-1) * x_i(2*i) - 50000, 2);
  return f;
}


