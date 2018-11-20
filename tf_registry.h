#ifndef TF_REGISTER
#define TF_REGISTER
#include <map>
#include <string>
#include "types.h"


using std::map;
using std::string;

struct FuncX0 {
  FuncX0(){}
  FuncX0(Function f, const vector<double>& p, int type) : _f(f), pattern(p), x0_type(type) {}
  Function _f;
  vector<double> pattern;
  int x0_type;
};

extern map<string,  FuncX0> tf_registry;

void init_registry();
#endif // TF_REGISTER

