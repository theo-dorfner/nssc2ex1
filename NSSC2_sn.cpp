#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <mpi.h>

/*
What I need from above:
  - NX, NY size of arrays (per process)
  - Anzahl der gesamt u 
  - solution, rhs, u
  - mean_runtime_it (mittlere runtime per Itteration)
  - anzahl der process 
  - falls es ein root rank gibt den sonst nehm ich 0
*/

double NormL2(const std::vector<long double> &v) {
  double norm = 0;
  for (const auto &value : v) {
    norm += value * value;
  }
  return sqrt(norm);
}

double NormInf(const std::vector<long double> &v) {
  double max = std::numeric_limits<double>::lowest();
  for (const auto &value : v) {
    max = std::fabs(value) > max ? std::fabs(value) : max;
  }
  return max;
}

int main(){

  int N = 4;
  int size = N*N;

  std::vector <long double> u(size, 1);                  // initialize solution arrays
  std::vector <long double> rhs(size, 3);                // intialise right hand side (size elements with initial value N)
  std::vector <double> solution(size, N);                // initialise solution

  double h = 1 / (1.0*N - 1);                            // nicht durch int dividieren --> N * 1.0
  double Center = 4 + h*h * 4*M_PI*M_PI;
  int NX = N, NY = N;

  // Initializations
  std::vector <long double> residual_elements(size, N);
  std::vector <long double> error_elemets(size, 0);
  double North, South, East, West;

  // Loop to calculate Residual and Error for NX*NY Matrix
  for(int i = 0; i < NX*NY; ++i){
    
    North = u[i-NX]; East  = u[i+1];
    South = u[i+NX]; West  = u[i-1];
  
    if(   i   % NX == 0     ){ West = 0.0;          // if no leftborder neighbour
    }if((i+1) % NX == 0     ){ East = 0.0;          // if no rightborder neighbour
    }if(  i   < NX          ){ North = 0.0;         // if no topborder neighbour
    }if(  i   > NX*(NY-1)-1 ){ South = 0.0;         // if no bottomborder neighbour
    }
    residual_elements[i] =  u[i]*Center - East - West - South - North - rhs[i]*h*h;
    error_elemets[i] = std::abs(u[i]-solution[i]);
  }

  // das hoff ich bekomme ich auch von oben 
  double mean_runtime_it; // Summe von allen Runtimes per Itteration divided by Itteration
  int process;  // Anzahl der thread

  // Combining all calculations to compute Residual und Error
  int root = 0
  if(myrank == root){

    // initializing buffer befor gather/sum
    std::vector <long double> residual(size of all u);
    std::vector <long double> error(size of all u);
    double runtime_sum;

    // gather all elements to combine
    MPI_Gather( &residual_elements.data(), residual_elements.size(), MPI_LONG_DOUBLE,
                residual.data(), residual_elements.size(), MPI_LONG_DOUBLE,
                root, MPI_COMM_WORLD);

    MPI_Gather( &error_elemets.data(), error_elemets.size(), MPI_LONG_DOUBLE,
                error.data(), error_elemets.size(), MPI_LONG_DOUBLE,
                root, MPI_COMM_WORLD);

    // sum up the mean-runtime and calculate mean
    MPI_Reduce( &mean_runtime_it, runtime_sum, 1, MPI_DOUBLE, 
                MPI_SUM, root, MPI_COMM_WORLD);
    double mean_runtime = runtime_sum/process;     
    
    // Calculating NormL2 and Infinite Norm of Residual and Error 
    auto residualNorm = NormL2(residual);
    auto residualMax = NormInf(residual);
    auto errorNorm = NormL2(error);
    auto errorMax = NormInf(error);
  
    // Output the result
    std::cout << std::scientific << "|residual|=" << residualNorm << std::endl;
    std::cout << std::scientific << "|residualMax|=" << residualMax << std::endl;
    std::cout << std::scientific << "|error|=" << errorNorm << std::endl;
    std::cout << std::scientific << "|errorMax|=" << errorMax << std::endl;
    std::cout << std::scientific << "average runtime per iteration = " << mean_runtime << std::endl;
  }

MPI_Finalize();
return 0;
  
}
