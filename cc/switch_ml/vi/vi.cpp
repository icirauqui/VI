#include "vi.h"







UGMM::UGMM(std::vector<double> X, int K, float sigma2) : _X(X), _N(X.size()), _K(K), _sigma2(sigma2) {
  std::cout << "UGMM constructor" << std::endl;
}


UGMM::~UGMM() {
  std::cout << "UGMM destructor" << std::endl;
}


void UGMM::Init() {
  std::cout << "UGMM Init" << std::endl;

  double scale = 10.0;

  std::mt19937 rng(0);
  std::normal_distribution<double> distribution(0.5, 0.5);

  double dirichlet_alpha = distribution(rng) * (int)(distribution(rng) * scale);

  std::vector<double> alpha(_N, dirichlet_alpha);
  _phi = std::vector<std::vector<double>>(_N, std::vector<double>(_K, 0.0));

  gsl_rng *r_global;
  gsl_rng_env_setup();
  r_global = gsl_rng_alloc(gsl_rng_default);
  for (int i = 0; i < _N; i++) {
    gsl_ran_dirichlet(r_global, _K, alpha.data(), _phi[i].data());
  }

  // Generate m
  double min = *std::min_element(_X.begin(), _X.end());
  double max = *std::max_element(_X.begin(), _X.end());

  // std::vector<int> m(_K);
  std::uniform_int_distribution<int> distribution_m(min, max);
  _m = std::vector<double>(_K, distribution_m(rng));
  for (int i = 0; i < _K; i++) {
    _m[i] += max * distribution(rng);
  }

  // Generate vector s2
  _s2 = std::vector<double>(_K, _sigma2 * distribution(rng));
}


float UGMM::GetELBO() {
  std::cout << std::endl << " - UGMM GetELBO" << std::endl;

  std::vector<double> t1(_K);
  for (int i = 0; i < _K; i++) {
    t1[i] = log(_s2[i]) - _m[i] / _sigma2;
  }

  double t1sum = std::accumulate(t1.begin(), t1.end(), 0.0);

  //print_vector("s2", _s2);
  //print_vector("m", _m);
  //std::cout << "sigma2: " << _sigma2 << std::endl;
  //print_vector("t1", t1);
  //std::cout << "t1sum: " << t1sum << std::endl;

  std::vector<double> outer_a(_N);
  for (int i = 0; i < _N; i++) {
    outer_a[i] = pow(_X[i], 2);
  }

  std::vector<double> outer_b(_K);
  for (int i = 0; i < _K; i++) {
    outer_b[i] = _s2[i] + pow(_m[i], 2);
  }

  std::vector<std::vector<double>> t2 = add_outer<double>(outer_a, outer_b, -0.5);

  std::vector<std::vector<double>> outer_X_m = product_outer<double>(_X, _m);

  for (unsigned int i = 0; i < t2.size(); i++) {
    for (unsigned int j = 0; j < t2[i].size(); j++) {
      t2[i][j] += outer_X_m[i][j];
      t2[i][j] -= log(_phi[i][j]);
      t2[i][j] *= _phi[i][j];
    }
  }

  double t2sum = 0.0;
  for (unsigned int i = 0; i < t2.size(); i++) {
    t2sum += std::accumulate(t2[i].begin(), t2[i].end(), 0.0);
  }

  //std::cout << "t2sum: " << t2sum << std::endl;

  std::cout << "      ELBO: " << t1sum + t2sum << std::endl;

  return t1sum + t2sum;
}



void UGMM::Fit(int max_iter, float tol) {
  std::cout << "UGMM Fit" << std::endl;

  //Init();

  std::vector<double> elbo_vals = std::vector<double>(1, GetELBO());
  std::vector<std::vector<double>> m_hist = std::vector<std::vector<double>>(1, _m);
  std::vector<std::vector<double>> s2_hist = std::vector<std::vector<double>>(1, _s2);

  for (int iter = 0; iter < max_iter; iter++) {

    // Print separator 
    std::cout << std::endl;
    std::cout << " - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
    std::cout << " - Iteration " << iter+1 << std::endl;
    std::cout << std::endl;

    Cavi();
    double elbo_new = GetELBO();

    elbo_vals.push_back(elbo_new);
    m_hist.push_back(_m);
    s2_hist.push_back(_s2);

    std::cout << std::endl << " - Result" << std::endl;

    if ((iter+1) % 5 == 0) {
        std::cout << "   iter: " << iter+1 << std::endl;
        std::cout << "      m: ";
        for (unsigned int j = 0; j<m_hist[iter].size(); j++) {
            std::cout << m_hist[iter][j] << " ";
        }
        std::cout << std::endl;
    }

    if (abs(elbo_new - elbo_vals[iter]) < tol) {
        std::cout << "      Converged with " << elbo_new << " at iteration " << iter << std::endl;
        break;
    }

    if (iter == max_iter - 1) {
        std::cout << "      Warning: max_iter reached, elbo " << elbo_new << std::endl;
    }

  }


}


void UGMM::Cavi() {
  UpdatePhi();
  UpdateMu();
}


void UGMM::UpdatePhi() {
  std::cout << std::endl << " - UGMM UpdatePhi" << std::endl;

  std::vector<std::vector<double>> t1 = product_outer<double>(_X, _m);

  std::vector<double> t2(_K);
  for (int i = 0; i < _K; i++) {
    t2[i] = -1 * (0.5 * pow(_m[i], 2) + 0.5 * _s2[i]);
  }

  std::vector<std::vector<double>> exponent = std::vector<std::vector<double>>(_N, std::vector<double>(_K, 0.0));
  for (int i = 0; i < _N; i++) {
    double sum = 0.0;
    for (int j = 0; j < _K; j++) {
      exponent[i][j] = t1[i][j] + t2[j];
      _phi[i][j] = std::exp(exponent[i][j]);
      sum += _phi[i][j];
    }
    if (sum==0.0)
      std::cout << "sum is 0" << std::endl;
    for (int j = 0; j < _K; j++) {
      _phi[i][j] /= sum;
    }
  }

  std::cout << "      phi: " << std::endl;
  for (int i = 0; i < 5; i++) {
    std::cout << "      ";
    for (int j = 0; j < _K; j++) {
      std::cout << _phi[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

void UGMM::UpdateMu() {
  std::cout << std::endl << " - UGMM UpdateMu" << std::endl;
  std::vector<std::vector<double>> phi_X = std::vector<std::vector<double>>(_N, std::vector<double>(_K, 0.0));
  std::vector<double> phi_X_sum = std::vector<double>(_K, 0.0);
  for (int i = 0; i < _N; i++) {
    for (int j = 0; j < _K; j++) {
      phi_X[i][j] = _phi[i][j] * _X[i];
      phi_X_sum[j] += phi_X[i][j];
    }
  }

  std::vector<double> phi_sum0_1 = std::vector<double>(_K, 0.0);
  for (int i = 0; i < _N; i++) {
    for (int j = 0; j < _K; j++) {
      phi_sum0_1[j] += _phi[i][j];
    }
  }

  for (int i = 0; i < _K; i++) {
    phi_sum0_1[i] = phi_sum0_1[i] + (1 / _sigma2);
    phi_sum0_1[i] = 1.0 / phi_sum0_1[i];
  }

  for (int i = 0; i < _K; i++) {
    _m[i] = phi_X_sum[i] * phi_sum0_1[i];
    _s2[i] = phi_sum0_1[i];
  }

  std::cout << "      m: " << std::endl << "      ";
  for (int i = 0; i < _K; i++) {
    std::cout << _m[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "      s2: " << std::endl << "      ";
  for (int i = 0; i < _K; i++) {
    std::cout << _s2[i] << " ";
  }
  std::cout << std::endl;

}



void UGMM::SetX(std::vector<double> X) {
  _X = X;
}

void UGMM::SetPhi(std::vector<std::vector<double>> phi) {
  _phi = phi;
}

void UGMM::SetM(std::vector<double> m) {
  _m = m;
}

void UGMM::SetS2(std::vector<double> s2) {
  _s2 = s2;
}

std::vector<double> UGMM::GetX() {
  return _X;
}

std::vector<std::vector<double>> UGMM::GetPhi() {
  return _phi;
}

std::vector<double> UGMM::GetM() {
  return _m;
}

std::vector<double> UGMM::GetS2() {
  return _s2;
}


