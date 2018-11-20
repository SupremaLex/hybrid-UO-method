#include "test_function.h"
#include "test_function_collection.h"
#include "tf_registry.h"


TestFunction::TestFunction(Function f, const int& n, const Vector& x0, const string& name) : _f(f), _dim(n), _nfev(0), _x0(x0), _name(name) {
  if (x0.size() != n)
    std::cout << "WTF" << std::endl;
}


double TestFunction::operator ()(const Vector& x) {
  _nfev++;
  return _f(x);
}


int TestFunction::get_dimension() const {
  return _dim;
}


int TestFunction::get_nfev() const {
  return _nfev;
}

Vector TestFunction::get_x0() const {
  return _x0;
}


string TestFunction::get_name() const {
  return _name;
}

TestFunctionFactory::TestFunctionFactory() {
  init_registry();
}

TestFunction TestFunctionFactory::create(const std::string& name,int dimension) {
  FuncX0 data  = tf_registry[name];
  return TestFunction(data._f, dimension, test_functions::create_x0(data.pattern, dimension, data.x0_type), name);
}

std::vector<std::string> TestFunctionFactory::get_all_names() {
  std::vector<std::string> names;
  names.reserve(tf_registry.size());
  for(auto const& pair: tf_registry)
      names.push_back(pair.first);
  return names;
}
