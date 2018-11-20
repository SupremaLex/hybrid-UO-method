#ifndef TEST_FUNCTION_COLLECTION
#define TEST_FUNCTION_COLLECTION
#include "types.h"
#define x_i(i) x(i-1)
#define sum(a, b, f) for (int i = a; i < b + 1; i++) f +=

namespace test_functions {

inline Vector create_x0(const vector<double>& pattern, int n, int type=0) {
  Vector x0 (n);
  int m = pattern.size();
  switch (type) {
  case 0:
    for (int i = 0; i < n; i++)
      x0(i) = pattern[i % m];
    break;
  case 1: {
    for (int i = 1; i < n + 1; i++)
      x0(i-1) = i;
  }
    break;
  case 2: {
    double h = 1 / (n + 1);
    for (int i = 1; i < n + 1; i++)
      x0(i-1) = i * h;
  }
    break;
  case 3: {
    for (int i = 1; i < n + 1; i++)
      x0(i-1) = 1.0 / i;
  }
    break;
  case 4:
    for (int i = 1; i < n + 1; i++)
      x0(i-1) = i / n;
    break;
  case 5:
    for (int i = 1; i < n + 1; i++)
      x0(i-1) = i / (n + 1);
    break;
  case 6: {
    for (int i = 1; i < n + 1; i++)
      x0(i-1) = 2.0;
    x0(0) = 1;
  }
    break;
  case 7:
    for (int i = 1; i < n + 1; i++)
      x0(i-1) = 1 - i / (n + 1);
    break;
  }
    return x0;
}


namespace extended {
  double ext_freudenstein_roth(const Vector& x);
  double ext_trigonometric(const Vector& x);
  double ext_rosenbrock(const Vector& x);
  double ext_white_holst(const Vector& x);
  double ext_beale(const Vector& x);
  double ext_penalty(const Vector& x);
  double ext_tet(const Vector& x);
  double ext_himmelblau(const Vector& x);
  double ext_PSC1(const Vector& x);
  double ext_powell(const Vector& x);
  double ext_maratos(const Vector& x);
  double ext_cliff(const Vector& x);
  double ext_wood(const Vector& x);
  double ext_hiebert(const Vector& x);
}  // namespace extended

namespace generalized {
  double gen_rosenbrock(const Vector& x);
  double gen_white_holst(const Vector& x);
  double gen_PSC1(const Vector& x);
}  // namespace generalized

namespace diagonal {
  double diagonal_1(const Vector& x);
  double diagonal_2(const Vector& x);
  double diagonal_3(const Vector& x);
  double diagonal_4(const Vector& x);
  double diagonal_5(const Vector& x);
  double diagonal_6(const Vector& x);
  double diagonal_7(const Vector& x);
  double diagonal_8(const Vector& x);
  double diagonal_9(const Vector& x);
  double ext_tridiagonal(const Vector& x);
  double ext_block_diagonal_1(const Vector& x);
  double ext_tridiagonal_2(const Vector& x);
  double broyden_tridiagonal(const Vector& x);
  double gen_tridiagonal_1(const Vector& x);
  double gen_tridiagonal_2(const Vector& x);
}  // namespace diagonal

namespace quadratic {
  double quadratic_1(const Vector& x);
  double quadratic_2(const Vector& x);
  double perturbed_quadratic(const Vector& x);
  double perturbed_quadratic_diagonal(const Vector& x);
  double ext_quadratic_penalty_1(const Vector& x);
  double ext_quadratic_penalty_2(const Vector& x);
  double ext_quadratic_exponential_1(const Vector& x);
  double partial_perturbed_quadratic(const Vector& x);
  double almost_perturbed_quadratic(const Vector& x);
  double perturbed_tridiagonal_quadratic(const Vector& x);
  double gen_quadratic(const Vector& x);
}  // namespace quadratic

namespace full_hessian {
  double full_hessian_1(const Vector& x);
  double full_hessian_2(const Vector& x);
  double full_hessian_3(const Vector& x);
}  // namespace full_hessian

namespace cute {
  double FLETCBV3(const Vector& x);
  double FLETCHCR(const Vector& x);
  double BDQRTIC(const Vector& x);
  double TRIDIA(const Vector& x);
  double ARWHEAD(const Vector& x);
  double NONDIA(const Vector& x);
  double NONDQUAR(const Vector& x);
  double DQDRTIC(const Vector& x);
  double EG2(const Vector& x);
  double LIARWHD_1(const Vector& x);
  double POWER(const Vector& x);
  double ENGVAL1(const Vector& x);
  double CRAGGLVY(const Vector& x);
  double EDENSCH(const Vector& x);
  double INDEF(const Vector& x);
  double CUBE(const Vector& x);
  double BDEXP(const Vector& x);
  double HARKERP2(const Vector& x);
  double GENHUMPS(const Vector& x);
  double MCCORMCK(const Vector& x);
  double NONSCOMP(const Vector& x);
  double VARDIM(const Vector& x);
  double QUARTC(const Vector& x);
  double SINQUAD(const Vector& x);
  double LIARWHD_2(const Vector& x);
  double DIXON3DQ(const Vector& x);
  double SINCOS(const Vector& x);
  double COSINE(const Vector& x);
  double SINE(const Vector& x);
  double BIGGSB1(const Vector& x);
  double HIMMELBG(const Vector& x);
  double HIMMELH(const Vector& x);
  double ext_DENSCHNB(const Vector& x);
  double ext_DENSCHNF(const Vector& x);
}  // namespace cube
}  // namespace test_functions

#endif // TEST_FUNCTION_COLLECTION

