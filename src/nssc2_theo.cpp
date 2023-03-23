#include <vector>
#include "jacobi.hpp"

using matrixType = std::vector <std::vector <double>>;

/*
    things needed
    matrix A
    initional rhs (boundary) b
    dimensions NX, NY
    number of iterations
*/

int main(){
    //variable declaration
    int fullSize = NX * NY;
    double solutionU[2][fullSize];
    double ghostValues[fullSize];
    double ghostInNorth[NX], ghostInSouth[NX]; // are dimensional allocations correct here?
    double ghostOutNorth[NX], ghostOutSouth[NX];
    int counter;
    int idNorth, idSouth;
    // prepare runtime

    // start iterations
    for(counter = 0; counter < iterations; ++counter){
        // need to collect runtime

        // collect neighbour IDs
        //int MPI_Cart_shift(MPI_Comm comm_cart, int direction, int disp, int *rank_source, int *rank_dest)
        int MPI_Cart_shift(MPI_Comm comm_cart, 1, -1, int *rank_source, int *rank_dest) //vertical - north
        int MPI_Cart_shift(MPI_Comm comm_cart, 1, +1, int *rank_source, int *rank_dest) //vertical - south

        // initiate non-blocking receive
        //MPI_Irecv( buf, count, datatype, source, tag, comm, [OUT] &request_handle);
        MPI_Irecv( &ghostInNorth, NY, MPI_Double, source, counter, comm, [OUT] &request_handle);
        MPI_Irecv( &ghostInSouth, NY, MPI_Double, source, tag, comm, [OUT] &request_handle);

        // initiate send

        // wait for receive

        // --- start jacobi-calculation
        // resolve neighbouring arrays into fullSize index position

        // actually calculate jacobi
    }

    // prepare output
    double finalSolution[fullSize];
    double rhs[fullSize];
    for(int i=0; i < fullSize; ++i){
        finalSolution[i] = solutionU[counter % 2][i];
        rhs[i] = b[i] + ghostValues[i];
    }
    // calculate mean runtime

}