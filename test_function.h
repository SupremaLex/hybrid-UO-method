#ifndef TEST_FUNCTION_H
#define TEST_FUNCTION_H
#include <vector>
#include <string>
#include "tf_registry.h"
#include "types.h"


struct TestFunction {
  TestFunction(Function f, const int& n, const Vector& x0, const string& name);
  double operator()(const Vector& x);
  int get_dimension() const;
  int get_nfev() const;
  Vector get_x0() const;
  string get_name() const;

 private:
  Function _f;
  int _dim;
  int _nfev;
  Vector _x0;
  string _name;
};

struct TestFunctionFactory {
  TestFunctionFactory();
  TestFunction create(const std::string& name, int dimension);
  std::vector<std::string> get_all_names();
};

#endif // TEST_FUNCTION_H
