#ifndef HYBRID_METHOD
#define HYBRID_METHOD
#include <math.h>
#include <vector>
#include <deque>
#include <utility>
#include <iostream>
#include "types.h"
#include "test_function.h"


struct TestResult {
    TestResult(const int& nit_, const int& nfev_, const Vector& x_, const double& min_, const double& grad_norm_) :
      nit(nit_), nfev(nfev_), x(x_), min(min_), grad_norm(grad_norm_) {}
    int nit;
    int nfev;
    Vector x;
    double min;
    double grad_norm;
};

/**
 * hybrid_method represents hybrid of quasi-Newton methods and p steps method
 */
TestResult hybrid_method(TestFunction f, int p, int method, double eps, int max_iter);
/**
 * gradient calculate function gradient by finite differenct (central)
 */
Vector gradient(TestFunction tf, Vector& x0, double eps);

void update_H(const deque<Vector>& x_k, const deque<Vector>& s_k,
              const Matrix& I, Matrix* h, int method = 0);
/**
 * DFP Davidon–Fletcher–Powell formula for Hessian update
 */
void DFP(const deque<Vector>& x_k, const deque<Vector>& s_k, const Matrix& I, Matrix* h);
/**
 *BFGS Broyden–Fletcher–Goldfarb–Shanno formula for Hessian update
 */
void BFGS(const deque<Vector>& x_k, const deque<Vector>& s_k, const Matrix& I, Matrix* h);
/**
 * calculate_beta inexact line search whith strong Wolfe condition
 */
std::pair<double, double> calculate_beta(TestFunction f, const Vector& x,
                                         const Vector& s, const Vector& grad, double eps);
Vector calc_gammas(const deque<Vector>& grads, int p, const Matrix& h);
bool stop_criteria(TestFunction f, const Vector& x1, const Vector& x2, const Vector& grad, double eps);
bool any_nan(const Matrix& m);
#endif // HYBRID_METHOD

