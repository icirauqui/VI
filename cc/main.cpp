#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <sstream>

#include "cc/switch_ml/vi/vi_agent.h"
#include "cc/switch_ml/vi/vi.h"

#include <chrono>

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


std::vector<std::vector<double>> ReadCSV(std::string path, char delimiter = ';') {
  std::ifstream file(path);
  std::vector<std::vector<double>> dataList;

  std::string line = "";
  // Iterate through each line and split the content using delimeter
  while (getline(file, line)) {
    std::vector<double> vec;
    std::stringstream ss(line);
    std::string cell;

    while (getline(ss, cell, delimiter)) {
      std::replace( cell.begin(), cell.end(), ',', '.'); 
      vec.push_back(std::stod(cell));
    }
    dataList.push_back(vec);
  }

  return dataList;
}



std::vector<std::vector<double>> vector_float_to_double(std::vector<std::vector<float>> X) {
  std::vector<std::vector<double>> result;
  for (unsigned int i=0; i<X.size(); i++) {
    std::vector<double> Xi;
    for (unsigned int j=0; j<X[i].size(); j++) {
      Xi.push_back(X[i][j]);
    }
    result.push_back(Xi);
  }
  return result;
}

std::vector<double> vector_2d_to_1d(std::vector<std::vector<double>> X) {
  std::vector<double> result;
  for (unsigned int i=0; i<X.size(); i++) {
    for (unsigned int j=0; j<X[i].size(); j++) {
      result.push_back(X[i][j]);
    }
  }
  return result;
}


void print_vec_1d(std::string name, std::vector<std::vector<double>> v) {
  std::cout << std::endl;
  std::cout << name << " [" << v.size() << "][" << v[0].size() << "]" << std::endl;
  for (unsigned int i=0; i<v.size(); i++) {
    for (unsigned int j=0; j<v[i].size(); j++) {
      std::cout << "  " << v[i][j];
    }
  }
  std::cout << std::endl;
}

void print_vec_2d(std::string name, std::vector<std::vector<double>> v, unsigned int limit = 5) {
  std::cout << std::endl;
  std::cout << name << " [" << v.size() << "][" << v[0].size() << "]" << std::endl;
  for (unsigned int i=0; i<limit; i++) {
    for (unsigned int j=0; j<v[i].size(); j++) {
      std::cout << "  " << v[i][j];
    }
    std::cout << std::endl;
  }
}


void print_separator() {
  std::cout << std::endl;
  std::cout << "--------------------------------------------------" << std::endl;
  std::cout << std::endl;
}

int main() {
    
  // Read data
  std::vector<std::vector<double>> m   = ReadCSV("/home/icirauqui/w4rkspace/VI/data/m.csv");
  std::vector<std::vector<double>> phi = ReadCSV("/home/icirauqui/w4rkspace/VI/data/phi.csv");
  std::vector<std::vector<double>> s2  = ReadCSV("/home/icirauqui/w4rkspace/VI/data/s2.csv");
  std::vector<std::vector<double>> X   = ReadCSV("/home/icirauqui/w4rkspace/VI/data/X.csv");

  print_vec_1d("m", m);
  print_vec_2d("phi", phi);
  print_vec_1d("s2", s2);
  print_vec_2d("X", X, 5);


  // Convert to double
  //std::vector<std::vector<double>> m_double   = vector_float_to_double(load_m);
  //std::vector<std::vector<double>> phi_double = vector_float_to_double(phi);
  //std::vector<std::vector<double>> s2_double  = vector_float_to_double(load_s2);
  //std::vector<std::vector<double>> X_double   = vector_float_to_double(load_X);

  std::vector<double> m_1d = vector_2d_to_1d(m);
  //std::vector<float> phi_1d = vector_2d_to_1d(m);
  std::vector<double> X_1d = vector_2d_to_1d(X);
  std::vector<double> s2_1d = vector_2d_to_1d(s2);

  print_separator();

  auto start1 = std::chrono::high_resolution_clock::now();

  int K = 3;
  UGMM ugmm(X_1d, K, 1.0);
  ugmm.Init();

  // Set simulation values from Python
  ugmm.SetPhi(phi);
  ugmm.SetS2(s2_1d);
  ugmm.SetM(m_1d);

  
  /*

  //Get vectors
  std::vector<std::vector<double>> get_phi  = ugmm.GetPhi();
  std::vector<double> get_m                 = ugmm.GetM();
  std::vector<double> get_s2                = ugmm.GetS2();

  //Compare put and get values
  std::cout << std::endl;
  std::cout << "m" << std::endl;
  for (unsigned int i=0; i<m_1d.size(); i++) {
    std::cout << "  " << m_1d[i] << "/" << get_m[i];
  }
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "s2" << std::endl;
  for (unsigned int i=0; i<s2_1d.size(); i++) {
    std::cout << "  " << s2_1d[i] << "/" << get_s2[i];
  }
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "phi" << std::endl;
  for (unsigned int i=0; i<5; i++) {
    for (unsigned int j=0; j<phi[i].size(); j++) {
      std::cout << "  " << phi[i][j] << "/" << get_phi[i][j];
    }
    std::cout << std::endl;
  }

  */
  


  print_separator();



  auto start12 = std::chrono::high_resolution_clock::now();

  ugmm.Fit();

  auto stop1 = std::chrono::high_resolution_clock::now();
  auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1);

  auto duration12 = std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start12);

  double elbo = ugmm.GetELBO();
  std::cout << std::endl << std::endl << "Final ELBO = " << elbo << std::endl;

  print_separator();

  auto start2 = std::chrono::high_resolution_clock::now();
  vi::Fit(phi, s2_1d, m_1d,
          X_1d, K, 1.0);
  auto stop2 = std::chrono::high_resolution_clock::now();
  auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);

  print_separator();

  std::cout << std::endl << std::endl;
  std::cout << "Time 1   = " << duration1.count() << " microseconds" << std::endl;
  std::cout << "Time 1.2 = " << duration12.count() << " microseconds" << std::endl;
  std::cout << "Time 2   = " << duration2.count() << " microseconds" << std::endl;


  //ugmm.Fit();
  //std::cout << std::endl << std::endl << "Init" << std::endl;
  //ugmm.Init();
  //std::cout << std::endl << std::endl << "GetELBO" << std::endl;
  //ugmm.GetELBO();









  //ugmm.Init();
  //ugmm.GetELBO();

  
  //for (unsigned int i=0; i<X.size(); i++) {
  //  std::cout << X[i] << " ";
  //}
  //std::cout << std::endl;

  //hello_vi();

  return 0;





}





/*

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





  // Read data
  std::vector<std::vector<float>> load_m   = ReadCSV("/home/icirauqui/w4rkspace/VI/data/m.csv");
  std::vector<std::vector<float>> load_phi = ReadCSV("/home/icirauqui/w4rkspace/VI/data/phi.csv");
  std::vector<std::vector<float>> load_s2  = ReadCSV("/home/icirauqui/w4rkspace/VI/data/s2.csv");
  std::vector<std::vector<float>> load_X   = ReadCSV("/home/icirauqui/w4rkspace/VI/data/X.csv");

  std::cout << "m size   " << load_m.size() << " " << load_m[0].size() << std::endl;
  std::cout << "phi size " << load_phi.size() << " " << load_phi[0].size() << std::endl;
  std::cout << "s2 size  " << load_s2.size() << " " << load_s2[0].size() << std::endl;
  std::cout << "X size   " << load_X.size() << " " << load_X[0].size() << std::endl;

  std::cout << "m" << std::endl;
  for (unsigned int i=0; i<load_m.size(); i++) {
    for (unsigned int j=0; j<load_m[i].size(); j++) {
      std::cout << "  " << load_m[i][j];
    }
  }
  std::cout << std::endl;

  std::cout << "s2" << std::endl;
  for (unsigned int i=0; i<load_s2.size(); i++) {
    for (unsigned int j=0; j<load_s2[i].size(); j++) {
      std::cout << "  " << load_s2[i][j];
    }
  }
  std::cout << std::endl;

  std::cout << "phi" << std::endl;
  for (unsigned int i=0; i<5; i++) {
    for (unsigned int j=0; j<load_phi[i].size(); j++) {
      std::cout << "  " << load_phi[i][j];
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;



  // Convert to double
  //std::vector<std::vector<double>> m   = vector_float_to_double(load_m);
  std::vector<std::vector<double>> phi = vector_float_to_double(load_phi);
  //std::vector<std::vector<double>> s2  = vector_float_to_double(load_s2);
  //std::vector<std::vector<double>> X   = vector_float_to_double(load_X);

  std::vector<float> X_1d = vector_2d_to_1d(load_X);
  std::vector<float> m_1d = vector_2d_to_1d(load_m);
  std::vector<float> s2_1d = vector_2d_to_1d(load_s2);

  UGMM ugmm(X_1d, K, 1.0);
  ugmm.SetPhi(phi);
  ugmm.SetS2(s2_1d);
  ugmm.SetM(m_1d);

  //Get vectors
  std::vector<std::vector<double>> get_phi = ugmm.GetPhi();
  std::vector<float> get_m    = ugmm.GetM();
  std::vector<float> get_s2   = ugmm.GetS2();

  //Compare put and get values
  std::cout << std::endl;
  std::cout << "m" << std::endl;
  for (unsigned int i=0; i<m_1d.size(); i++) {
    std::cout << "  " << m_1d[i] << "/" << get_m[i];
  }
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "s2" << std::endl;
  for (unsigned int i=0; i<s2_1d.size(); i++) {
    std::cout << "  " << s2_1d[i] << "/" << get_s2[i];
  }
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "phi" << std::endl;
  for (unsigned int i=0; i<5; i++) {
    for (unsigned int j=0; j<phi[i].size(); j++) {
      std::cout << "  " << phi[i][j] << "/" << get_phi[i][j];
    }
    std::cout << std::endl;
  }




  //ugmm.Fit();
  //std::cout << std::endl << std::endl << "Init" << std::endl;
  //ugmm.Init();
  //std::cout << std::endl << std::endl << "GetELBO" << std::endl;
  //ugmm.GetELBO();









  //ugmm.Init();
  //ugmm.GetELBO();

  
  //for (unsigned int i=0; i<X.size(); i++) {
  //  std::cout << X[i] << " ";
  //}
  //std::cout << std::endl;

  //hello_vi();

  return 0;
}

*/