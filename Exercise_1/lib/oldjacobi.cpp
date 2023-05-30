#include <iostream>
#include <vector>

using matrixType = std::vector <std::vector <double>>;

std::vector <double> jacobi(matrixType A,std::vector <double> b,std::vector <double> x,int dim){
    std::vector <double> nX(dim);
    for(int i=0; i<dim; i++){
        double sum{0};
        for(int j=0; j<dim; j++) {
            if(i==j) continue;
            sum += A[i][j]*x[j];
        }
        nX[i] = (b[i] - sum)/A[i][i];
    }
    return nX;
}
void printMatrix(matrixType matrix,int N){
    for(int i=0; i < N; i++){
        for(int j=0; j < N; j++){
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}
void printVector(std::vector <double> vector){
    for(auto& elem:vector) std::cout << elem << " ";
    std::cout << std::endl;
}
/*
int main(){
    int nN = 2;
    std::vector <std::vector <double>> ol(nN,std::vector <double> (nN,2));
    std::vector <std::vector <double>> a{{2,1},{-1,2}};
    std::vector <double> b{1,3};
    std::vector <double> x(1,0);
    printMatrix(a,nN);
    for(int i=0; i<10; i++){
        x = jacobi(a,b,x,nN);
        printVector(x);
    }
    return 0;
}
*/

