#include "test_function_collection.h"
namespace cute = test_functions::cute;
#include "test_function.h"

double cute::FLETCBV3(const Vector& x) {
  /*Another Boundary Value problem*/
  int n = x.size();
  double p = 1e-8;
  double f = 0.5 * p * (pow(x_i(1), 2) + pow(x_i(n), 2));
  double h = 1 / (n + 1);
  double hh = h * h;
  sum(1, n-1, f)
      0.5 * p * pow(x_i(i) - x_i(i+1), 2) -
      x_i(i) * p * (hh + 2) / hh + cos(x_i(i)) * p / hh;
  return f - x_i(n) * p * (hh + 2) / hh + cos(x_i(n)) * p / hh;
}


double cute::FLETCHCR(const Vector& x) {
    /*Another Boundary Value problem*/
  double f = 0;
  int n = x.size();
  sum(1, n - 1, f)
      100 * pow(x_i(i+1) - x_i(i) + 1 - pow(x_i(i), 2), 2);
  return f;
}


double cute::BDQRTIC(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n - 4, f)
      pow(-4 * x_i(i) + 3, 2) +
      pow(pow(x_i(i), 2) + 2 * pow(x_i(i+1), 2) + 3 * pow(x_i(i+2), 2) + 4 * pow(x_i(i+3), 2) + 5 * pow(x_i(n), 2) , 2);
  return f;
}


double cute::TRIDIA(const Vector& x) {
  double f = pow(x_i(1) - 1, 2);
  int n = x.size();
  sum(2, n, f)
      pow(2 * x_i(i) - x_i(i-1), 2) * i;
  return f;
}


double cute::ARWHEAD(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n - 1, f)
      (-4 * x_i(i) + 3) +
      pow(pow(x_i(i), 2) + pow(x_i(n), 2), 2);
  return f;
}


double cute::NONDIA(const Vector& x) {
  double f = pow(x_i(1) - 1, 2);
  int n = x.size();
  sum(2, n, f)
      100 * pow(x_i(i) - pow(x_i(i-1), 2), 2);
  return f;
}


double cute::NONDQUAR(const Vector& x) {
  int n = x.size();
  double f = pow(x_i(1) - x_i(2), 2) + pow(x_i(n-1) - x_i(n), 2);
  sum(1, n - 2, f)
      pow(x_i(i) + x_i(i+1) + x_i(n), 4);
  return f;
}


double cute::DQDRTIC(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n - 2, f)
      pow(x_i(i), 2) + 100 * pow(x_i(i+1), 2) + 100 *pow(x_i(i+2), 2);
  return f;
}


double cute::EG2(const Vector& x) {
  int n = x.size();
  double f = pow(sin(x_i(n)), 2);
  sum(1, n - 1, f)
      sin(x_i(1) + pow(x_i(i), 2) - 1);
  return f;
}


double cute::LIARWHD_1(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n, f)
      4 * pow(-x_i(1) + pow(x_i(i), 2), 2) + pow(x_i(i) - 1, 2);
  return f;
}


double cute::POWER(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n, f)
      pow(i * x_i(i), 2);
  return f;
}


double cute::ENGVAL1(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n - 1, f)
      pow(pow(x_i(i), 2) + pow(x_i(i+1), 2), 2) +
      (-4 * x_i(i) + 3);
  return f;
}


double cute::CRAGGLVY(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n / 2 - 1, f)
      pow(exp(x_i(2*i-1)) - x_i(2*i), 4) + 100 * pow(x_i(2*i) - x_i(2*i+1), 6) +
      pow(tan(x_i(2*i+1) + x_i(2*i + 2)) + x_i(2*i+1) - x_i(2*i + 2), 4) + pow(x_i(2*i-1), 8) +
      pow(x_i(2*i+2) - 1, 2);
  return f;
}


double cute::EDENSCH(const Vector& x) {
  double f = 16;
  int n = x.size();
  sum(1, n - 1, f)
      pow(x_i(i) - 2, 4) + pow(x_i(i) * x_i(i+1) - 2 * x_i(i+1), 2) +
      pow(x_i(i+1) + 1, 2);
  return f;
}


double cute::INDEF(const Vector& x) {
  /* indefinite Hessian at start point*/
  double f = 0;
  int n = x.size();
  sum(1, n - 1, f)
      x_i(i) + 0.5 * cos(2 * x_i(i) - x_i(n) - x_i(1));
  return f + x_i(n);
}


double cute::CUBE(const Vector& x) {
  double f = pow(x_i(1) - 1, 2);
  int n = x.size();
  sum(2, n, f)
      100 * pow(x_i(i) - pow(x_i(i-1), 3), 2);
  return f;
}


double cute::BDEXP(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n - 2, f)
      (x_i(i) + x_i(i+1)) * exp(-x_i(i+2) * (x_i(i) + x_i(i+1)));
  return f;
}

double cute::HARKERP2(const Vector& x) {
  double f1 = 0, f2 = 0, f3 = 0;
  int n = x.size();
  sum(1, n, f1)
      x_i(i);
  sum(1, n, f2)
      x_i(i) + 0.5 * pow(x_i(i), 2);
  for(int i = 2; i < n + 1; i++) {
    double tmp = 0;
    for (int j = i; j < n + 1; j++)
      tmp +=x_i(j);
    f3 += tmp * tmp;
  }
  return f1 * f1 - f2 + 2 * f3;
}


double cute::GENHUMPS(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n - 1, f)
      pow(sin(2 * x_i(i)), 2) * pow(sin(2 * x_i(i+1)), 2) +
      0.05 * (pow(x_i(i), 2) + pow(x_i(i+1), 2));
  return f;
}


double cute::MCCORMCK(const Vector& x) {
  /* -> - inf*/
  double f = 0;
  int n = x.size();
  sum(1, n - 1, f)
      -1.5 * x_i(i) + 2.5 * x_i(i+1) + 1 + pow(x_i(i) - x_i(i+1), 2) +
      sin(x_i(i) + x_i(i+1));
  return f;
}


double cute::NONSCOMP(const Vector& x) {
  double f = pow(x_i(1) - 1, 2);
  int n = x.size();
  sum(2, n, f)
      4 * pow(x_i(i) - pow(x_i(i-1), 2), 2);
  return f;
}


double cute::VARDIM(const Vector& x) {
  double f1 = 0, f2 = 0;
  int n = x.size();
  sum(1, n, f1)
      pow(x_i(i) - 1, 2);
  sum(1, n, f2)
      i * x_i(i) - n * (n + 1) / 2;
  return f1 + f2 * f2 + f2 * f2 * f2 *f2;
}


double cute::QUARTC(const Vector& x) {
  int n = x.size();
  double f = 0;
  sum(1, n, f)
      pow(x_i(i) - 1, 4);
  return f;
}


double cute::SINQUAD(const Vector& x) {
  int n = x.size();
  double f = pow(x_i(1) - 1, 4) + pow(pow(x_i(n), 2) - pow(x_i(1), 2), 2);
  sum(2, n - 1, f)
      pow(sin(x_i(i) - x_i(n)) - pow(x_i(1), 2) + pow(x_i(i), 2), 2);
  return f;
}


double cute::LIARWHD_2(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n, f)
      4 * pow(pow(x_i(i), 2) - x_i(1), 2) + pow(x_i(i) - 1, 2);
  return f;
}


double cute::DIXON3DQ(const Vector& x) {
  int n = x.size();
  double f = pow(x_i(1) - 1, 2) + pow(x_i(n) - 1, 2);
  sum(1, n - 1, f)
      pow(x_i(i) - x_i(i+1), 2);
  return f;
}


double cute::SINCOS(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n / 2, f)
      pow(pow(x_i(2*i-1), 2) + pow(x_i(2*i), 2) + x_i(2*i-1) * x_i(2*i), 2) +
      pow(sin(x_i(2*i-1)), 2) + pow(cos(x_i(2*i)), 2);
  return f;
}


double cute::COSINE(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n - 1, f)
      cos(-0.5 * x_i(i+1) + pow(x_i(i), 2));
  return f;
}


double cute::SINE(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n - 1, f)
      sin(-0.5 * x_i(i+1) + pow(x_i(i), 2));
  return f;
}


double cute::BIGGSB1(const Vector& x) {
  int n = x.size();
  double f = pow(x_i(1) - 1, 2) + pow(1 - x_i(n), 2);
  sum(1, n - 1, f)
      pow(x_i(i+1) - x_i(i), 2);
  return f;
}


double cute::HIMMELBG(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n / 2, f)
      (2 * pow(x_i(2*i-1), 2) + 3 * pow(x_i(2*i), 2)) * exp(-x_i(2*i-1) - x_i(2*i));
  return f;
}


double cute::HIMMELH(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n / 2, f)
      (-3 * x_i(2*i-1) - 2 * x_i(2*i) + 2 + pow(x_i(2*i-1), 3) + pow(x_i(2*i), 2));
  return f;
}


double cute::ext_DENSCHNB(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n / 2, f)
      pow(x_i(2*i-1) - 2, 2) * (1 + pow(x_i(2*i), 2)) + pow(x_i(2*i) + 1, 2);
  return f;
}


double cute::ext_DENSCHNF(const Vector& x) {
  double f = 0;
  int n = x.size();
  sum(1, n / 2, f)
      pow(2 * pow(x_i(2*i-1) + x_i(2*i), 2) + pow(x_i(2*i-1) - x_i(2*i), 2) - 8, 2) +
      pow(5 * pow(x_i(2*i-1), 2) + pow(x_i(2*i) - 3, 2) - 9, 2);
  return f;
}
