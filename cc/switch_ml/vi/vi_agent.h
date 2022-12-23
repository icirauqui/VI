#ifndef VI_AGENT_H
#define VI_AGENT_H

#include <iostream>
#include <vector>
#include <random>
#include <string>
#include <sstream>
#include <fstream>

#include <algorithm>
#include <cstdlib>
#include <cmath>

// GLS Dirichlet
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>




namespace vi {

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



std::vector<std::vector<double>> InitPhi(std::vector<double> &X, int N, int K);
std::vector<double> InitM(std::vector<double> &X, int N, int K, float sigma2);
std::vector<double> InitS2(std::vector<double> &X, int N, int K, float sigma2);

double GetELBO(std::vector<std::vector<double>> &phi, std::vector<double> &m, std::vector<double> &X, std::vector<double> &s2, int N, int K, float sigma2);

void Fit(std::vector<double> &_X, int _K, float _sigma2, int max_iter = 100, float tol = 1e-10);

void Fit(std::vector<std::vector<double>> &phi, std::vector<double> &s2, std::vector<double> &m,
         std::vector<double> &_X, int _K, float _sigma2, int max_iter = 100, float tol = 1e-10);

void UpdatePhi(std::vector<std::vector<double>> &phi, std::vector<double> &m, std::vector<double> &X, std::vector<double> &s2, int N, int K, float sigma2);

void UpdateMu(std::vector<std::vector<double>> &phi, std::vector<double> &m, std::vector<double> &X, std::vector<double> &s2, int N, int K, float sigma2);




} // namespace vi

#endif // VI_AGENT_H