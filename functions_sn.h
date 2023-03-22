#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>


std::vector <long double> Residual_Calc( int NX, int NY,
    std::vector <long double> &u, std::vector <long double> &rhs)
{

    std::vector <long double> residual_elements(NX*NY, 0);      // Initializations
    double North, South, East, West;
    double h = 1 / (1.0*NX - 1);                            // nicht durch int dividieren --> N * 1.0
    double Center = 4 + h*h * 4*M_PI*M_PI;

    for(int i = 0; i < NX*NY; ++i)                    // Loop to calculate Residual
    {
      North = u[i-NX]; East  = u[i+1];                // Set values
      South = u[i+NX]; West  = u[i-1];
    
      if(   i   % NX == 0     ){ West = 0.0;          // if no leftborder neighbour
      }if((i+1) % NX == 0     ){ East = 0.0;          // if no rightborder neighbour
      }if(  i   < NX          ){ North = 0.0;         // if no topborder neighbour
      }if(  i   > NX*(NY-1)-1 ){ South = 0.0;         // if no bottomborder neighbour
      }
      residual_elements[i] =  u[i]*Center - East - West - South - North - rhs[i]*h*h;
    }
    return residual_elements;

}


std::vector <long double> Error_Calc(int NX, int NY,
    std::vector <long double> &u, std::vector <long double> &solution)
{   
    std::vector <long double> error_elemets(NX*NY, 0);  // Initialization

    for(int i = 0; i < NX*NY; ++i)
    {error_elemets[i] = std::abs(u[i]-solution[i]);     // Loop to calculate error
    }
    return error_elemets;

}


double NormL2(const std::vector<long double> &v) 
{

  double norm = 0;
  for (const auto &value : v) {
    norm += value * value;
  }
  return sqrt(norm);

}


double NormInf(const std::vector<long double> &v) 
{

  double max = std::numeric_limits<double>::lowest();
  for (const auto &value : v) {
    max = std::fabs(value) > max ? std::fabs(value) : max;
  }
  return max;

}