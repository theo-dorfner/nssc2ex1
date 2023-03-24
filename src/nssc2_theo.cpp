#include <vector>
#include "jacobi.hpp"
#include <chrono>
#include <cmath>
# include <mpi.h>

using matrixType = std::vector <std::vector <double>>;

using Clock = std::chrono::system_clock;

/*
    things needed:
    matrix A
    initional rhs (boundary) b
    dimensions NX, NY
    number of iterations
*/

int main(){
    //variable declaration
    int fullSize = NX * NY;
    double solutionU[2][fullSize]{0};
    double ghostValues[fullSize];
    double ghostInNorth[NX], ghostInSouth[NX]; // are dimensional allocations correct here?
    double ghostOutNorth[NX], ghostOutSouth[NX];
    int idNorth, idSouth, procID;
    std::chrono::duration<double> procRuntime{0};
    MPI_Request requestNorth;
    MPI_Request requestSouth;
    MPI_Status statusNorth;
    MPI_Status statusSouth;


    // collect neighbour IDs
    //int MPI_Cart_shift(MPI_Comm comm_cart, int direction, int disp, int *rank_source, int *rank_dest)
    MPI_Cart_shift(MPI_Comm comm_cart, 1, -1, &procID, &idNorth); //vertical - north
    MPI_Cart_shift(MPI_Comm comm_cart, 1, +1, &procID, &idSouth); //vertical - south

    // start iterations
    for(int counter = 0; counter < iterations; ++counter){
        auto start = std::chrono::steady_clock::now(); // start runtime timing

        // generate outgoing ghost layers
        for(int i=0; i<NX;++i) ghostOutNorth[i] = solutionU[(counter-1)%2][(NY-1)*NX+i];
        for(int i=0; i<NX;++i) ghostOutSouth[i] = solutionU[(counter-1)%2][i];

        // initiate non-blocking receive
        //MPI_Irecv( buf, count, datatype, source, tag, comm, [OUT] &request_handle);
        MPI_Irecv( &ghostInNorth, NY, MPI_Double, idNorth, counter, comm, &requestNorth);
        MPI_Irecv( &ghostInSouth, NY, MPI_Double, idSouth, counter, comm, &requestSouth);

        // initiate send
        MPI_Send(&ghostOutNorth, NX,MPI_Double, idNorth, counter, MPI_Comm_cart);
        MPI_Send(&ghostOutSouth, NX,MPI_Double, idSouth, counter, MPI_Comm_cart);

        // wait for receive
        MPI_Wait(&requestNorth,&statusNorth);
        MPI_Wait(&requestSouth,&statusSouth);

        // --- start jacobi-calculation
        // resolve neighbouring arrays into fullSize index position
        for(int i=0; i<NX;++i) g[(NY-1)*NX+i] = (+1) * ghostInNorth[i]; //need to clarify +1/-1 here
        for(int i=0; i<NX;++i) g[i] = (+1) * ghostInSouth[i];

        // actually calculate jacobi
        for(int i=0; i<dim; i++){
            double sum{0};
            for(int j=0; j<fullSize; j++) {
                if(i==j) continue;
                sum += A[i][j] * solutionU[(counter+1)%2][i];
            }
            solutionU[counter%2][i] = (b[i] + g[i] - sum)/A[i][i];
        }

        //calc runtime
        procRuntime += std::chrono::steady_clock::now() - start;
    }

    //if(counter == iterations) std::cout << "problem with counter not being interations" << std::endl;

    // prepare output
    double finalSolution[fullSize];
    double rhs[fullSize];
    for(int i=0; i < fullSize; ++i){
        finalSolution[i] = solutionU[iterations % 2][i];
        rhs[i] = b[i] + ghostValues[i];
    }
    // calculate mean runtime in seconds
    double meanRuntime = procRuntime.count()/(std::pow(10,9)*iterations*1.0);

}// A[j*N+i] - i = NX; j = NY