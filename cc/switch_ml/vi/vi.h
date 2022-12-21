#ifndef VI_H
#define VI_H

#include <iostream>
#include <vector>
#include <random>
#include <string>

// Dirichlet distribution
//#include <gls/gls_rng.h>
#include <gsl/gsl_linalg.h>
// GLS Dirichlet
//#include <gls/gls_dirichlet.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


//#include <gsl/gsl_math.h>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_test.h>

//#include <gsl/gsl_ieee_utils.h>

//#include <gsl/gsl_integration.h>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_blas.h>
//#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_cdf.h>
//#include <gsl/gsl_statistics.h>

#include <algorithm>
#include <cstdlib>
#include <cmath>

#include <sstream>
#include <fstream>

template <typename T>
void print_vector(std::string name, std::vector<T> v) {
  std::cout << name << ": ";
  for (unsigned int i = 0; i < v.size(); i++) {
      std::cout << v[i] << " ";
  }
  std::cout << std::endl;
}


template <typename T1, typename T2, typename T3>
std::vector<std::vector<T1>> product_outer(std::vector<T2> &a, std::vector<T3> &b, float factor = 1.0) {
  std::vector<std::vector<T1>> c(a.size(), std::vector<T1>(b.size(), 0.0));
  for (unsigned int i = 0; i < a.size(); i++) {
      for (unsigned int j = 0; j < b.size(); j++) {
          c[i][j] = a[i] * b[j];
      }
  }
  return c;
}


template <typename T1, typename T2, typename T3>
std::vector<std::vector<T1>> add_outer(std::vector<T2> &a, std::vector<T3> &b, float factor = 1.0) {
  std::vector<std::vector<T1>> c(a.size(), std::vector<T1>(b.size(), 0.0));
  for (unsigned int i = 0; i < a.size(); i++) {
      for (unsigned int j = 0; j < b.size(); j++) {
          c[i][j] = factor * (a[i] + b[j]);
      }
  }
  return c;
}








class UGMM {

public:

  UGMM(std::vector<float> X, int K, float sigma2);

  ~UGMM();

  void Init();

  float GetELBO();

  void Fit(int max_iter = 100, float tol = 1e-10);

  void Cavi();

  void UpdatePhi();

  void UpdateMu();





  void SetPhi(std::vector<std::vector<double>> phi);

  void SetM(std::vector<float> m);
  
  void SetS2(std::vector<float> s2);

  std::vector<std::vector<double>> GetPhi();

  std::vector<float> GetM();

  std::vector<float> GetS2();




private:

  std::vector<float> _X;

  std::vector<std::vector<double>> _phi;
  
  int _N;
  int _K;
  float _sigma2;

  std::vector<float> _m;
  std::vector<float> _s2;





};



#endif // VI_H