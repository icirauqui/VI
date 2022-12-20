#include "vi.h"



UGMM::UGMM(std::vector<float> X, int K, float sigma2):
    _X(X), _N(X.size()), _K(K), _sigma2(sigma2) {
    std::cout << "UGMM constructor" << std::endl;
}


UGMM::~UGMM() {
    std::cout << "UGMM destructor" << std::endl;
}


void UGMM::Init() {
    float scale = 10.0;

    std::mt19937 rng(0);
    std::normal_distribution<float> distribution(0.5, 0.5);

    float dirichlet_alpha = distribution(rng) * (int)(distribution(rng) * scale);

    std::vector<double> alpha(_N, dirichlet_alpha);
    std::vector<std::vector<double>> phi(_N, std::vector<double>(_K, 0.0));

    gsl_rng *r_global;
    gsl_rng_env_setup ();
    r_global = gsl_rng_alloc (gsl_rng_default);
    for (int i=0; i<_N; i++) {
        gsl_ran_dirichlet (r_global, _K, alpha.data(), phi[i].data());
    }

    // Generate m
    float min = *std::min_element(_X.begin(), _X.end());
    float max = *std::max_element(_X.begin(), _X.end());
    
    //std::vector<int> m(_K);
    std::uniform_int_distribution<int> distribution_m(min, max);
    _m = std::vector<int>(_K, distribution_m(rng));
    for (int i = 0; i < _K; i++) {
        _m[i] += max * distribution(rng);
    }

    // Generate vector s2
    _s2 = std::vector<float>(_K, _sigma2 * distribution(rng));
}


float UGMM::GetELBO() {
    std::cout << "UGMM GetELBO" << std::endl;

    std::vector<float> t1(_K);
    for (int i = 0; i < _K; i++) {
        t1[i] = log(_s2[i]) - _m[i] / _sigma2;
    }
    
    float t1sum = std::accumulate(t1.begin(), t1.end(), 0.0);

    print_vector("s2", _s2);
    print_vector("m", _m);
    std::cout << "sigma2: " << _sigma2 << std::endl;
    print_vector("t1", t1);


    return 0.0;
}


void UGMM::Fit() {
    std::cout << "UGMM Fit" << std::endl;
    Init();


}