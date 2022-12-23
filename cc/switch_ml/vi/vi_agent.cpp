#include "vi_agent.h"







std::vector<std::vector<double>> vi::InitPhi(std::vector<double> &X, int N, int K) {
  std::cout << "UGMM Init Phi" << std::endl;

  double scale = 10.0;

  std::mt19937 rng(0);
  std::normal_distribution<double> distribution(0.5, 0.5);

  double dirichlet_alpha = distribution(rng) * (int)(distribution(rng) * scale);

  std::vector<double> alpha(N, dirichlet_alpha);
  std::vector<std::vector<double>> _phi = std::vector<std::vector<double>>(N, std::vector<double>(K, 0.0));

  gsl_rng *r_global;
  gsl_rng_env_setup();
  r_global = gsl_rng_alloc(gsl_rng_default);
  for (int i = 0; i < N; i++) {
    gsl_ran_dirichlet(r_global, K, alpha.data(), _phi[i].data());
  }

  return _phi;
}




std::vector<double> vi::InitM(std::vector<double> &X, int N, int K, float sigma2) {
  std::cout << "UGMM Init m" << std::endl;

  std::mt19937 rng(0);
  std::normal_distribution<double> distribution(0.5, 0.5);

  // Generate m
  double min = *std::min_element(X.begin(), X.end());
  double max = *std::max_element(X.begin(), X.end());

  // std::vector<int> m(_K);
  std::uniform_int_distribution<int> distribution_m(min, max);
  std::vector<double> _m = std::vector<double>(K, distribution_m(rng));
  for (int i = 0; i < K; i++) {
    _m[i] += max * distribution(rng);
  }

  return _m;
}




std::vector<double> vi::InitS2(std::vector<double> &X, int N, int K, float sigma2) {
  std::cout << "UGMM Init s2" << std::endl;

  std::mt19937 rng(0);
  std::normal_distribution<double> distribution(0.5, 0.5);

  std::vector<double> _s2 = std::vector<double>(K, sigma2 * distribution(rng));

  return _s2;
}





double vi::GetELBO(std::vector<std::vector<double>> &phi, std::vector<double> &m, std::vector<double> &X, std::vector<double> &s2, int N, int K, float sigma2) {
  std::cout << std::endl << " - UGMM GetELBO" << std::endl;

  std::vector<double> t1(K);
  for (int i = 0; i < K; i++) {
    t1[i] = log(s2[i]) - m[i] / sigma2;
  }

  double t1sum = std::accumulate(t1.begin(), t1.end(), 0.0);

  std::vector<double> outer_a(N);
  for (int i = 0; i < N; i++) {
    outer_a[i] = pow(X[i], 2);
  }

  std::vector<double> outer_b(K);
  for (int i = 0; i < K; i++) {
    outer_b[i] = s2[i] + pow(m[i], 2);
  }

  std::vector<std::vector<double>> t2 = add_outer<double>(outer_a, outer_b, -0.5);

  std::vector<std::vector<double>> outer_X_m = product_outer<double>(X, m);

  for (unsigned int i = 0; i < t2.size(); i++) {
    for (unsigned int j = 0; j < t2[i].size(); j++) {
      t2[i][j] += outer_X_m[i][j];
      t2[i][j] -= log(phi[i][j]);
      t2[i][j] *= phi[i][j];
    }
  }

  double t2sum = 0.0;
  for (unsigned int i = 0; i < t2.size(); i++) {
    t2sum += std::accumulate(t2[i].begin(), t2[i].end(), 0.0);
  }

  std::cout << "      ELBO: " << t1sum + t2sum << std::endl;

  return t1sum + t2sum;
}



void vi::Fit(std::vector<double> &_X, int _K, float _sigma2, int max_iter, float tol) {

  std::cout << "UGMM Fit" << std::endl;

  int _N = _X.size();

  std::vector<std::vector<double>> phi = vi::InitPhi(_X, _N, _K);
  std::vector<double> m = vi::InitM(_X, _N, _K, _sigma2);
  std::vector<double> s2 = vi::InitS2(_X, _N, _K, _sigma2);

  double elbo = GetELBO(phi, m, _X, s2, _N, _K, _sigma2);

  std::vector<double> hist_elbo = std::vector<double>(1, elbo);
  std::vector<std::vector<double>> hist_m = std::vector<std::vector<double>>(1, m);
  std::vector<std::vector<double>> hist_s2 = std::vector<std::vector<double>>(1, s2);

  for (int iter = 0; iter < max_iter; iter++) {
    // Print separator 
    std::cout << std::endl;
    std::cout << " - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
    std::cout << " - Iteration " << iter+1 << std::endl;
    std::cout << std::endl;

    // CAVI
    UpdatePhi(phi, m, _X, s2, _N, _K, _sigma2);
    UpdateMu(phi, m, _X, s2, _N, _K, _sigma2);

    double elbo_new = GetELBO(phi, m, _X, s2, _N, _K, _sigma2);

    hist_elbo.push_back(elbo_new);
    hist_m.push_back(m);
    hist_s2.push_back(s2);

    std::cout << std::endl << " - Result" << std::endl;

    if ((iter+1) % 5 == 0) {
        std::cout << "   iter: " << iter+1 << std::endl;
        std::cout << "      m: ";
        for (unsigned int j = 0; j<hist_m[iter].size(); j++) {
            std::cout << hist_m[iter][j] << " ";
        }
        std::cout << std::endl;
    }

    if (abs(elbo_new - hist_elbo[iter]) < tol) {
        std::cout << "      Converged with " << elbo_new << " at iteration " << iter << std::endl;
        break;
    }

    if (iter == max_iter - 1) {
        std::cout << "      Warning: max_iter reached, elbo " << elbo_new << std::endl;
    }

  }
}



void vi::Fit(std::vector<std::vector<double>> &phi, std::vector<double> &s2, std::vector<double> &m,
             std::vector<double> &_X, int _K, float _sigma2, int max_iter, float tol) {

  std::cout << "UGMM Fit" << std::endl;

  int _N = _X.size();

  double elbo = GetELBO(phi, m, _X, s2, _N, _K, _sigma2);

  std::vector<double> hist_elbo = std::vector<double>(1, elbo);
  std::vector<std::vector<double>> hist_m = std::vector<std::vector<double>>(1, m);
  std::vector<std::vector<double>> hist_s2 = std::vector<std::vector<double>>(1, s2);

  for (int iter = 0; iter < max_iter; iter++) {
    // Print separator 
    std::cout << std::endl;
    std::cout << " - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
    std::cout << " - Iteration " << iter+1 << std::endl;
    std::cout << std::endl;

    // CAVI
    UpdatePhi(phi, m, _X, s2, _N, _K, _sigma2);
    UpdateMu(phi, m, _X, s2, _N, _K, _sigma2);

    double elbo_new = GetELBO(phi, m, _X, s2, _N, _K, _sigma2);

    hist_elbo.push_back(elbo_new);
    hist_m.push_back(m);
    hist_s2.push_back(s2);

    std::cout << std::endl << " - Result" << std::endl;

    if ((iter+1) % 5 == 0) {
        std::cout << "   iter: " << iter+1 << std::endl;
        std::cout << "      m: ";
        for (unsigned int j = 0; j<hist_m[iter].size(); j++) {
            std::cout << hist_m[iter][j] << " ";
        }
        std::cout << std::endl;
    }

    if (abs(elbo_new - hist_elbo[iter]) < tol) {
        std::cout << "      Converged with " << elbo_new << " at iteration " << iter << std::endl;
        break;
    }

    if (iter == max_iter - 1) {
        std::cout << "      Warning: max_iter reached, elbo " << elbo_new << std::endl;
    }

  }
}



void vi::UpdatePhi(std::vector<std::vector<double>> &phi, std::vector<double> &m, std::vector<double> &X, std::vector<double> &s2, int N, int K, float sigma2) {
  std::cout << std::endl << " - UGMM UpdatePhi" << std::endl;

  std::vector<std::vector<double>> t1 = product_outer<double>(X, m);

  std::vector<double> t2(K);
  for (int i = 0; i < K; i++) {
    t2[i] = -1 * (0.5 * pow(m[i], 2) + 0.5 * s2[i]);
  }

  std::vector<std::vector<double>> exponent = std::vector<std::vector<double>>(N, std::vector<double>(K, 0.0));
  for (int i = 0; i < N; i++) {
    double sum = 0.0;
    for (int j = 0; j < K; j++) {
      exponent[i][j] = t1[i][j] + t2[j];
      phi[i][j] = std::exp(exponent[i][j]);
      sum += phi[i][j];
    }
    if (sum==0.0)
      std::cout << "sum is 0" << std::endl;
    for (int j = 0; j < K; j++) {
      phi[i][j] /= sum;
    }
  }

  std::cout << "      phi: " << std::endl;
  for (int i = 0; i < 5; i++) {
    std::cout << "      ";
    for (int j = 0; j < K; j++) {
      std::cout << phi[i][j] << " ";
    }
    std::cout << std::endl;
  }
}



void vi::UpdateMu(std::vector<std::vector<double>> &phi, std::vector<double> &m, std::vector<double> &X, std::vector<double> &s2, int N, int K, float sigma2) {
  std::cout << std::endl << " - UGMM UpdateMu" << std::endl;
  std::vector<std::vector<double>> phi_X = std::vector<std::vector<double>>(N, std::vector<double>(K, 0.0));
  std::vector<double> phi_X_sum = std::vector<double>(K, 0.0);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      phi_X[i][j] = phi[i][j] * X[i];
      phi_X_sum[j] += phi_X[i][j];
    }
  }

  std::vector<double> phi_sum0_1 = std::vector<double>(K, 0.0);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      phi_sum0_1[j] += phi[i][j];
    }
  }

  for (int i = 0; i < K; i++) {
    phi_sum0_1[i] = phi_sum0_1[i] + (1 / sigma2);
    phi_sum0_1[i] = 1.0 / phi_sum0_1[i];
  }

  for (int i = 0; i < K; i++) {
    m[i] = phi_X_sum[i] * phi_sum0_1[i];
    s2[i] = phi_sum0_1[i];
  }

  std::cout << "      m: " << std::endl << "      ";
  for (int i = 0; i < K; i++) {
    std::cout << m[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "      s2: " << std::endl << "      ";
  for (int i = 0; i < K; i++) {
    std::cout << s2[i] << " ";
  }
  std::cout << std::endl;

}




