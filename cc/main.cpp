#include <iostream>
#include <vector>
#include <random>
#include <fstream>

#include "cc/switch_ml/vi/vi_agent.h"
#include "cc/switch_ml/vi/vi.h"

template <typename T>
std::vector<T> slicing(std::vector<T>& arr, int X, int Y) {
 
    // Starting and Ending iterators
    auto start = arr.begin() + X;
    auto end = arr.begin() + Y + 1;
 
    // To store the sliced vector
    std::vector<T> result(Y - X + 1);
 
    // Copy vector using copy function()
    std::copy(start, end, result.begin());
 
    // Return the final sliced vector
    return result;
}

std::vector<float> normal_distribution(int N, float mean, float std, float scale = 1.0) {
  std::mt19937 rng(0);
  std::normal_distribution<float> distribution(mean, std);

  std::vector<float> result;
  for (int i = 0; i < N; i++) {
      result.push_back(scale * distribution(rng));
  }

  return result;
}


int write_to_file(std::vector<float> X, std::string filename) {
  //std::ofstream file;
  //file.open(filename);
  //for (unsigned int i=0; i<X.size(); i++) {
  //  file << X[i] << " ";
  //}
  //file.close();

  // Open the file
  std::FILE* fp = std::fopen("data.txt", "w");
  if (!fp) {
    std::cerr << "Error: Could not open file\n";
    return 1;
  }

  // Write the data to the file
  for (double x : X) {
    std::fprintf(fp, "%lf\n", x);
  }

  // Close the file
  std::fclose(fp);

  return 0;
}


int main() {
    
  int K = 3;
  int N = 1000;

  std::vector<float> mu = normal_distribution(K, 0.0, 1.0, 10.0);

  std::vector<float> X;
  for (unsigned int i=0; i<mu.size(); i++) {
    std::vector<float> Xi = normal_distribution(N, mu[i], 1.0, 10.0);
    X.insert(X.end(), Xi.begin(), Xi.end());
  }

  write_to_file(X, "/home/icirauqui/w4rkspace/VI/cc/switch_ml/vi/X.txt");

  // Plot data using gnuplot
  //std::system("gnuplot -p -e 'set term png; set output \"/home/icirauqui/w4rkspace/VI/cc/switch_ml/vi/X.png\"; set style data histograms; set style fill solid; plot \"/home/icirauqui/w4rkspace/VI/cc/switch_ml/vi/X.txt\" using 3'");
  
  std::vector<std::vector<float>> Xs;
  for (unsigned int i=0; i<mu.size(); i++) {
    std::vector<float> Xi = slicing(X, i*N, (i+1)*N-1);
    Xs.push_back(Xi);
  }


  std::ofstream file;
  std::string filename = "/home/icirauqui/w4rkspace/VI/cc/switch_ml/vi/X.txt";
  file.open(filename);
  file << "X1;X2;X3" << std::endl;
  for (unsigned int i=0; i<Xs[0].size(); i++) {
    for (unsigned int j=0; j<Xs.size(); j++) {
      file << Xs[j][i] << ";";
    }
    file << std::endl;
  }
  file.close();

  UGMM ugmm(X, K, 1.0);
  ugmm.Fit();

  //ugmm.Init();
  //ugmm.GetELBO();

  
  //for (unsigned int i=0; i<X.size(); i++) {
  //  std::cout << X[i] << " ";
  //}
  //std::cout << std::endl;

  //hello_vi();

  return 0;
}

