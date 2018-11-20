#include "hybrid_method.h"
#include <numeric>
#include <functional>
#include <algorithm>

using arma::norm;
using arma::dot;


TestResult hybrid_method(TestFunction f, int p, int method, double eps, int max_iter) {
  deque<Vector> x_k, g_k, s_k;
  Vector x0 = f.get_x0();
  Vector g = gradient(f, x0, eps);
  int k = 1;
  int n = f.get_dimension();
  double fx_prev = f(x0);
  double fx_last = fx_prev;
  double eps_2 = sqrt(eps);
  double eps_3 = pow(eps, 1.0 / 3);
  Matrix I = arma::eye<Matrix>(n, n);
  Matrix H = I;

  x_k.push_back(x0);
  g_k.push_back(g);
  s_k.push_back(-g);

  std::pair<double, double> pair_beta1_f = calculate_beta(f, x0, s_k.back(), g_k.back(), eps);
  x_k.push_back(x0 +  pair_beta1_f.first * s_k.back());
  fx_last = pair_beta1_f.second;

  while (/*stop_criteria(f, x_k.end()[-2], x_k.end()[-1], g_k.back(), eps)*/
         fabs(fx_last - fx_prev) >= eps * (1 + fabs(fx_last)) ||
         norm(x_k.end()[-1] - x_k.end()[-2], 2) >= eps_2 * (1 + norm(x_k.back()), 2) ||
         norm(g_k.back(), 2) > eps_3 * (1 + fabs(fx_last))
         && k < max_iter) {
    k++;
    if (s_k.size() > p) {
      x_k.pop_front();
      s_k.pop_front();
      g_k.pop_front();
    }
    g_k.push_back(gradient(f, x_k.back(), eps));
    if (k % n)
      update_H(x_k, g_k, I, &H, method);
    else
      H = I;
    Vector s = H * (-g_k.back());
    Vector gammas = calc_gammas(g_k, p, H);
    for (int i = 0; i < gammas.size(); i++)
      s += gammas[i] * s_k.end()[-i-1];
    s_k.push_back(s);
    std::pair<double, double> pair_beta_f = calculate_beta(f, x_k.back(), s_k.back(), g_k.back(), eps);
    x_k.push_back(x_k.back() + pair_beta_f.first * s_k.back());
    fx_prev = fx_last;
    fx_last = pair_beta_f.second;
    //std::cout << f(x_k.back()) << std::endl;
  }
  return TestResult(k, f.get_nfev(), x_k.back(), f(x_k.back()), norm(g_k.back(), 2));
}


Vector gradient(TestFunction f, Vector& x0, double eps) {
  int n = x0.size();
  Vector grad (n);
  double e = pow(std::numeric_limits<double>::epsilon(), 1.0 / 3);
  double h = e;
  for (int i = 0; i < n; i++) {
    /* Dennis J. E. , Robert B. Schnabel 5.4-5.6*/
    if (fabs(x0(i)) > eps)
      h = fabs(e * x0(i));
    double tmp = x0(i) + h;
    h = tmp - x0(i);
    x0(i) += h;
    double f_p = f(x0);
    x0(i) -= 2 * h;
    double f_m = f(x0);
    grad(i) = (f_p - f_m) / (2 * h);
    x0(i) += h;
    grad(i) = (f_p - f_m) / (2 * h);
  }
  return grad;
}

void update_H(const deque<Vector>& x_k, const deque<Vector>& s_k,
              const Matrix& I, Matrix* h, int method) {
  using arma::trans;
  Vector s = x_k.end()[-1] - x_k.end()[-2];
  Vector y = s_k.end()[-1] - s_k.end()[-2];
  if (method == 0) {
    double tmp = dot(y, s);
    if (tmp) {
      *h -= *h * y * trans(y) * (*h) / dot(trans(y) * (*h), y) + s * trans(s) / tmp;
    }
    else
      *h = I;
  } else {
    double tmp = dot(y, s);
    if (tmp) {
      *h = (I - s * trans(y) / tmp) * (*h) * (I - y * trans(s) / tmp) + s * trans(s) / tmp;
    }
    else
      *h = I;
  }
}

void DFP(const deque<Vector>& x_k, const deque<Vector>& s_k, const Matrix& I, Matrix* h) {
  using arma::trans;
  Vector s = x_k.end()[-1] - x_k.end()[-2];
  Vector y = s_k.end()[-1] - s_k.end()[-2];
  double tmp = dot(y, s);
  if (tmp) {
    *h -= *h * y * trans(y) * (*h) / dot(trans(y) * (*h), y) + s * trans(s) / tmp;
  }
  else
    *h = I;
}


void BFGS(const deque<Vector>& x_k, const deque<Vector>& s_k, const Matrix& I, Matrix* h) {
  Vector s = x_k.end()[-1] - x_k.end()[-2];
  Vector y = s_k.end()[-1] - s_k.end()[-2];
  double tmp = dot(y, s);
  if (tmp) {
    *h = (I - s * trans(y) / tmp) * (*h) * (I - y * trans(s) / tmp) + s * trans(s) / tmp;
  }
  else
    *h = I;
}


double quadratic_interpolation_min(double f0, double f1, double g, double a, double b) {
  auto quad_int = [f0, f1, g](double z) -> double {return g * z + f0 + (f1 - f0 - g) * z * z;};
  double c = - 0.5 * g / (f1 - f0 - g);
  double q0 = 0, q1, q2;
  if (c < a || c > b)
    c = a;
  q0 = quad_int(c);
  q1 = quad_int(a);
  q2 = quad_int(b);
  double min =  std::min({q0, q1, q2});
  if (q0 == min)
      return c;
  else if (q1 == min)
    return a;
  else
    return b;
}


double cubic_interpolation_min(double f0, double f1, double g0, double g1, double a, double b) {
  double eta = 3 * (f1 - f0) - 2 * g0 - g1;
  double ksi = g0 + g1 - 2 * (f1 - f0);
  auto cubic_int = [f0, g0, eta, ksi](double z) -> double {return f0 + g0 * z + eta * z * z + ksi * z * z * z;};
  double d = sqrt(eta * eta - 3 * ksi * g0);
  double c1 = (-eta - d) / (3 * ksi);
  double c2 = (-eta + d) / (3 * ksi);
  double q0 = 0, q1 = 0, q2, q3;
  if (c1 < a || c1 > b)
    c1 = a;
  if (c2 < a || c2 > b)
    c2 = a;
  q0 = cubic_int(c1);
  q1 = cubic_int(c2);
  q2 = cubic_int(a);
  q3 = cubic_int(b);

  double min =  std::min({q0, q1, q2, q3});
  if (q0 == min)
      return c1;
  else if (q1 == min)
    return c2;
  else if (q2 == min)
    return a;
  else
    return b;
}


std::pair<double, double> calculate_beta(TestFunction f, const Vector& x,
                                         const Vector& s, const Vector& grad, double eps) {
  /* R. Fletcher Practical methods of optimization (p. 35, alg 2.6.4)
   * inexact line search with strong Wolfe condition and quadratic interpolation
   * */
  /* Wolfe condition constants*/
  double c1 = 0.1, c2 = 0.3;
  /* Bracketing algorithm factors*/
  double t2 = 0.1, t3 = 0.5;
  double a = 0, b = 1;
  double f0 = f(x);
  double g0 = dot(grad, s);
  double alpha, f_alpha;
  while (fabs(a - b) > eps) {
    double a_ = a + t2 * (b - a);
    double b_ = b - t3 * (b - a);
    Vector t = x + a_ * s;
    //Vector t1 = x + b_ * s;
    double g1 = dot(gradient(f, t, eps), s);
    //double g2 = dot(gradient(f, t1, eps), s);
    double f1 = f(x + a_ * x);
    double f2 = f(x + b_ * s);
    /* Choose alpha as a minimum of quadratic interpolation*/
    alpha = a_ + quadratic_interpolation_min(f1, f2, g1, 0, 1) * (b_ - a_);
    Vector x_alpha = x + alpha * s;
    f_alpha = f(x_alpha);
    /* Armijo rule*/
    if (f_alpha > f0 + c1 * alpha * g0 || f_alpha >= f(x + a * s)) {
      b = alpha;
    } else {
      Vector g = gradient(f, x_alpha, eps);
      double tmp = dot(g, s);
      if( (a - alpha) * tmp <= eps)
        return std::make_pair(alpha, f_alpha);
      /* Wolfe rule*/
      if (fabs(dot(g, s)) <= c2 * fabs(tmp))
        return std::make_pair(alpha, f_alpha);
      a = alpha;
      if ((b - a) * tmp >= 0)
        b = a;
    }
  }
  return std::make_pair(alpha, f_alpha);
}

Vector calc_gammas(const deque<Vector>& g_k, int p, const Matrix& h) {
  int n = g_k.size();
  int k = n-1;
  int t = n < p ? n : p;
  Vector gammas (t-1);
  for (int i = 1; i < t; i++) {
    gammas[i-1] = dot(h * g_k[k], g_k[k - i + 1] - g_k[k - i]) / pow(norm(g_k[k - i], 2), 2);
  }
  return gammas;
}


/*bool stop_criteria(TestFunction f, const Vector& x1, const Vector& x2, const Vector& grad, double eps) {
  double y1 = f(x1);
  double y2 = f(x2);
  double abs_y2 = fabs(y2);
  bool cond1 = fabs(y1 - y2) < eps * (1 + abs_y2);
  bool cond2 = norm_2(x1 - x2) < sqrt(eps) * (1 + norm_2(x2));
  bool cond3 = norm_2(grad) <= pow(eps, 1.0 / 3) * (1 + abs_y2);
  return !(cond1 && cond2 && cond3);

}*/



/*bool any_nan(const Matrix& m) {
  for (int i = 0; i < m.size1(); i++) {
    for (int j = 0; j < m.size2(); j++) {
      if (std::isnan(m(i, j)))
        return true;
    }
  }
    return false;
}*/


