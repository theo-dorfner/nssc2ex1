#include <iostream>
#include "cppJacobi.hpp"
//#include <numbers>
#include <cmath>
#include <chrono>
#include <vector>

using Clock = std::chrono::system_clock;
using Duration = Clock::duration;

using std::cout;
using std::endl;
using std::vector;

double k = 2*M_PI;                // initialize k

// allow for easy function input
long double funcF(double x, double y){
    return (k*k*std::sin(x*2*M_PI)*std::sinh(y*2*M_PI)) ;
    //return 0;
}

int main(int argc, char **args){
    // retrieve program input
    if(argc != 3){
        cout << "wrong number of arguments - program terminated" << endl;
        return 99;
    }
    int N = strtol(args[1],nullptr,0);                                                  //resolution is N
    int iterations = strtol(args[2],nullptr,0);

    
    double dh = 1.0/(N - 1);                                                            // calculate spacing
    int size = N*N;
    vector <vector <long double>> u(2,vector<long double> (size,0));                    // initialize solution arrays
    vector <long double> rhs(size,0);                                                   // intialise right hand side
    vector <double> solution(size,0);                                                   // initialise solution

    for(int m=0;m < size;m++){
        double x{(m%N)*dh},y{(m/N)*dh};
        solution[m] = std::sin(k*x)*std::sinh(k*y);                                                         // construct solution
        if((1 - m/N/(N-1)) * (1 - m%N/(N-1)) * (1 - (N - m/(N))/(N)) * (1 - (N - m%(N))/(N)) == 0){         // filter out boundary values
            if(m/N/(N-1) == 1) {                                                                            // filter out "top-row-positions" (x,y=1)
                u[0][m] = std::sin(k*x) * std::sinh(k);
                u[1][m] = std::sin(k*x) * std::sinh(k);
            }
        }
        else{
        rhs[m] = funcF(x,y) * dh * dh;                                                                      // construct right hand side
        }
    }



    //perform iterations
    auto start = std::chrono::steady_clock::now();  //start time count

    for(int i = 0;i < iterations;i++){
        int toFill{(i+1)%2}, source{i%2};                                                                           // determine which array to use as source and which to write on
        for(int m=0;m < size;m++){
            if((1 - m/N/(N-1)) * (1 - m%N/(N-1)) * (1 - (N - m/(N))/(N)) * (1 - (N - m%(N))/(N)) == 0) continue;    // filter if current value is boundary value
            double temp = rhs[m];
            temp += u[source][m+N];
            temp += u[source][m+1];
            temp += u[source][m-1];
            temp += u[source][m-N];
            temp *= 1/(4.0 + dh*dh*k*k);
            u[toFill][m] = temp;                                // apply calculated value
        }
    }

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration = end - start;

    int final = iterations % 2;                                 // find vector from last iteration
    vector<double> error(size,0),residual(size,0);
    
    for(int m=0;m < size;m++){
        error[m] = std::abs(u[final][m] - solution[m]);
        if((1 - m/N/(N-1)) * (1 - m%N/(N-1)) * (1 - (N - m/(N))/(N)) * (1 - (N - m%(N))/(N)) == 0){     // for boundary values matrix entries are only diagonal elements
            residual[m] = u[final][m];
        }
        else{
            residual[m] = (4 + dh*dh*k*k)*u[final][m] - (u[final][m-N] + u[final][m-1] + u[final][m+1] + u[final][m+N]);
        }
        residual[m] = std::abs(residual[m] - rhs[m]); 
    }

    cout << "maximum total error: " << maximumNorm(error) << endl;
    cout << "euclidian norm of total error: " << euclidianNorm(error) << endl;
    cout << "maximum residual error: " << maximumNorm(residual) << endl;
    cout << "euclidian residual norm: " << euclidianNorm(residual) << endl;
    cout << "runtime: " << duration.count()/std::pow(10,9) << " s" << endl;

    return 0;
}