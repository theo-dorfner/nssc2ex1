#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>


std::vector <double> Residual_Calc( int NX, int NY, int resolution,
    std::vector <double> &finalSolution, std::vector <double> &rhs)
{

    std::vector <double> residual_elements(NX*NY, 0);      // Initializations
    double North, South, East, West;
    double h = 1.0 / (resolution - 1);                            // nicht durch int dividieren --> N * 1.0
    double Center = 4 + h*h * 4*M_PI*M_PI;

    for(int i = 0; i < NX*NY; ++i)                    // Loop to calculate Residual
    {
      North = finalSolution[i-NX]; East  = finalSolution[i+1];                // Set values
      South = finalSolution[i+NX]; West  = finalSolution[i-1];
    
      if(   i   % NX == 0     ){ West = 0.0;          // if no leftborder neighbour
      }if((i+1) % NX == 0     ){ East = 0.0;          // if no rightborder neighbour
      }if(  i   < NX          ){ North = 0.0;         // if no topborder neighbour
      }if(  i   > NX*(NY-1)-1 ){ South = 0.0;         // if no bottomborder neighbour
      }
      residual_elements[i] =  finalSolution[i]*Center - East - West - South - North - rhs[i];
    }
    return residual_elements;

}


std::vector <double> Error_Calc(int NX, int NY,
    std::vector <double> &u, std::vector <double> &solution)
{   
    std::vector <double> error_elemets(NX*NY, 0);  // Initialization

    for(int i = 0; i < NX*NY; ++i)
    {error_elemets[i] = std::abs(u[i]-solution[i]);     // Loop to calculate error
    }
    return error_elemets;

}

// L2 Norm ^ 2
double NormL2_2(const std::vector<double> &v)
{

  double sqnorm = 0;
  for (const auto &value : v) {
    sqnorm += value * value;
  }
  return sqnorm;

}


double NormInf(const std::vector<double> &v) 
{

  double max = std::numeric_limits<double>::lowest();
  for (const auto &value : v) {
    max = std::fabs(value) > max ? std::fabs(value) : max;
  }
  return max;

}