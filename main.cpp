#include <iostream>
#include <ctime>
#include <cstdlib>
#include <omp.h>
#include "cstdio"
#include <hybrid_method.h>
#include <test_function.h>
using namespace std;


int main(int argc, char** argv) {
  double eps = 1e-6;

  TestFunctionFactory factory = TestFunctionFactory();
  auto lst = factory.get_all_names();
  for (auto s : lst) {
      cout << s << endl;
      TestFunction f = factory.create(s, 100);
      auto res = hybrid_method(f, 2, 1, eps, 100000);
      cout << res.min << " "  << res.nit <<  "  " << res.nfev << endl;
  }


  vector<int> steps = {2};//, 3, 10, 50, 100, 200};
  vector<int> dims = {100};
  vector<int>::const_iterator it;
  clock_t begin = clock();
  double start = omp_get_wtime();
  int n = 1;
//#pragma omp parallel for
  for (it = steps.begin(); it < steps.end(); it++) {
    TestFunction f = factory.create("GeneralizedTridiagonal2", 1000);
    cout << f(f.get_x0()) << endl;
    /*for(int i = 0; i < n; i++) {
      auto res = hybrid_method(f, *it, 1, eps, 100000);
      cout << res.get_min() << " "  << res.get_nit() << endl;
      //cout << res.get_x() << endl;
    }*/
  }
  double time = omp_get_wtime() - start;
  cout << time / n << endl;
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout <<  elapsed_secs << endl;
  return 0;
}
/*Indef, Fletcv3, fletchcr,ExtendedHiebert ExtendedDENSCHNF EDENSCH*/
