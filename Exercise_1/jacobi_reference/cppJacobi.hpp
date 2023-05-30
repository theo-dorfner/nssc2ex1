#ifndef header_h
#define header_h

#include <vector>
#include <cmath>

using std::vector;
double maximumNorm(const vector <double> &vec){
    double result{0};
    for(auto &element:vec)if(element > result)result = element;
    return result;
}

double euclidianNorm(const vector <double> &vec){
    double result{0};
    for(auto &element:vec)result += element*element;
    return std::sqrt(result);
}


#endif