#include <mpi.h>
#include <stdio.h>
#include <iostream>

#include <cmath>
#include <vector>
#include <array>
#include <limits>
#include <chrono>

#include "../include/functions.h"
#include "../include/functions_sn.h"
#include "../include/jacobi.hpp"

using matrixType = std::vector <std::vector <double>>;

using Clock = std::chrono::system_clock;

int main(int argc, char* argv[]) {
    int ndims = 1;
    //int size; - had to comment out because it's redeclared later on
    int proc;
    int my_rank;
    MPI_Comm comm1D;
    int dims[ndims], coord;
    int wrap_around[ndims];
    int reorder;
    int my_cart_rank;
    int ierr;
    int nrows, ncols;

    MPI_Init(&argc, &argv);

    int resolution = atoi(argv[1]);
    int precs = resolution - 2;
    int iterations = atoi(argv[2]);


    MPI_Comm_size(MPI_COMM_WORLD, &proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    dims[0] = 0;

    MPI_Dims_create(proc, ndims, dims);

    wrap_around[0] = 0;
    reorder = 1;
    ierr = 0;
    ierr = MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, wrap_around, reorder, &comm1D);

    if(ierr != 0) printf("ERROR[%i] creating CART\n", ierr);

    //this is now done multiple times, but don't know how else to do it!
    std::vector<int>UPP;
    UPP = UnknownsPerProc(UPP, precs, proc); //calculates how big dimensions of arrays are for each rank!

    //IMPLEMENT FUNCTION(PROOF IF ENOUGH PRECS FOR PROCS)

    std::vector<double>b;
    std::vector<double>u(UPP[my_rank],0);
    std::vector<double>A(UPP[my_rank]*UPP[my_rank],0); //this is a Matrix --> initialize with A[j*N+i]
    //double* A[UPP[my_rank]]; //matrix; initialize with A[j*N + i]

    //u = Initialize_u0(u, UPP[my_rank]); //for "normal" application
    //b = Initialize_b0(b, UPP[my_rank]);
    double h = H(resolution);

    A = Initialize_A0(A, UPP[my_rank], precs, h);

    std::vector<int>y_begin;
    y_begin = UPPtoYBegin(UPP, precs, proc);

    std::vector<double>Y_vals;
    b = Initialize_b0(b, y_begin, precs, h, my_rank, proc);

    //vector_printer_int(y_begin);

    std::cout<< "rank: " << my_rank << std::endl;
    vector_printer(b);
    std::cout << std::endl;











// PART THEO
    //transition
    int &NX = ncols;
    int &NY = nrows;


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
    int startNorth = (NY-1)*NX;


    // collect neighbour IDs
    //int MPI_Cart_shift(MPI_Comm comm_cart, int direction, int disp, int *rank_source, int *rank_dest)
    MPI_Cart_shift(comm1D, 1, -1, &procID, &idNorth); //vertical - north
    MPI_Cart_shift(comm1D, 1, +1, &procID, &idSouth); //vertical - south

    // start iterations
    for(int counter = 0; counter < iterations; ++counter){
        auto start = std::chrono::steady_clock::now(); // start runtime timing

        // generate outgoing ghost layers
        for(int i=0; i<NX;++i) ghostOutNorth[i] = solutionU[(counter-1)%2][startNorth + i];
        for(int i=0; i<NX;++i) ghostOutSouth[i] = solutionU[(counter-1)%2][i];

        // initiate non-blocking receive
        //MPI_Irecv( buf, count, datatype, source, tag, comm, [OUT] &request_handle);
        MPI_Irecv( &ghostInNorth, NY, MPI_Double, idNorth, counter, comm1D, &requestNorth);
        MPI_Irecv( &ghostInSouth, NY, MPI_Double, idSouth, counter, comm1D, &requestSouth);

        // initiate send
        MPI_Send(&ghostOutNorth, NX,MPI_Double, idNorth, counter, comm1D);
        MPI_Send(&ghostOutSouth, NX,MPI_Double, idSouth, counter, comm1D);

        // wait for receive
        MPI_Wait(&requestNorth,&statusNorth);
        MPI_Wait(&requestSouth,&statusSouth);

        // --- start jacobi-calculation
        // resolve neighbouring arrays into fullSize index position
        for(int i=0; i<NX;++i) ghostValues[startNorth + i] = (+1) * ghostInNorth[i]; //need to clarify +1/-1 here
        for(int i=0; i<NX;++i) ghostValues[i] = (+1) * ghostInSouth[i];

        // actually calculate jacobi
        for(int i=0; i<fullSize; i++){
            double sum{0};
            for(int j=0; j<fullSize; j++) {
                if(i==j) continue;
                sum += A[i][j] * solutionU[(counter+1)%2][i];
            }
            solutionU[counter%2][i] = (b[i] + ghostValues[i] - sum)/A[i][i];
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


    /* things being passed on from my part:
    fullSize = NX*NY (this is the length of solution and rhs vector)

    finalSolution (of length fullSize) [afaik this is also the u that they want for everything else]
    rhs (this is b + ghostValues)
    meanRuntime




    */


// PART STEFFI

    int N = 4;
    int size = N*N; // sure that this shouldn't be NX*NY?

    std::vector <long double> u(size, 1);                  // initialize solution arrays
    std::vector <long double> rhs(size, 3);                // intialise right hand side (size elements with initial value N)
    std::vector <double> solution(size, N);                // initialise solution

    double h = 1 / (1.0*N - 1);                            // nicht durch int dividieren --> N * 1.0
    double Center = 4 + h*h * 4*M_PI*M_PI;
    int NX = N, NY = N;

    // Initializations
    std::vector <long double> residual_elements(NX*NY, 0);
    std::vector <long double> error_elemets(NX*NY, 0);
    
    // Calculate Residual and Error
    residual_elements = Residual_Calc(NX, NY, u, rhs);
    error_elemets = Error_Calc(NX, NY, u, solution);

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
