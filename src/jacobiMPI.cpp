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

    //A. MPI_Cart creation

    //A.1. Reading in values
    int resolution = atoi(argv[1]);
    int precs = resolution - 2;
    int iterations = atoi(argv[2]);

    //A.2. Reading out proc and rank; create cart and check if successful
    MPI_Comm_size(MPI_COMM_WORLD, &proc);
    //saves number of processes on "proc"
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    //saves current rank on "my_rank"

    dims[0] = 0;

    MPI_Dims_create(proc, ndims, dims);

    wrap_around[0] = 0;
    reorder = 1;
    ierr = 0;
    ierr = MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, wrap_around, reorder, &comm1D);
    //creating MPI_cart with one dimension

    if(ierr != 0) printf("ERROR[%i] creating CART\n", ierr);
    //returns error if cart creation wasn't successful

    //this is now done multiple times, but don't know how else to do it!
    std::vector<int>UPP;
    UPP = UnknownsPerProc(UPP, precs, proc); //calculates how big dimensions of arrays are for each rank!
    int NY = UPP[my_rank] / precs;
    //Check if resolution high enough for PROCS
    if(precs < proc)
    {
        if(my_rank == 0) cout << "Please choose a resolution that for N procs is at least N+2" << endl;
        return 0;
    }

    //B.1. Definition of Matrix A, vector u and vector b; size of arrays corresponds to UPP of the rank;
    vector<double>b;
    vector<double>u(UPP[my_rank],0);
    //starting value for vector u is all zero
    vector<double>A(UPP[my_rank]*UPP[my_rank],0);
    //this is a Matrix --> initialize with A[j*N+i]


    //B.2. Helpfunctions for the Initialization of A and b;
    double h = H(resolution);
    //calculating h
    vector<int>y_begin;
    //vector for the initialization of b
    y_begin = UPPtoYBegin(UPP, precs, proc);

    //B.3. Initialization of A and b;
    A = Initialize_A0(A, UPP[my_rank], precs, h);
    b = Initialize_b0(b, y_begin, precs, h, my_rank, proc);

    /*
    cout << "Hello, my rank is " << my_rank << " and the number of unknowns are " << UPP[my_rank] << endl << endl;
    cout << "This is my vector b: " << endl;
    vector_printer(b);
    cout <<"And this is my matrix A: " << endl;
    matrix_printer(A, UPP[my_rank]);
    cout<<"Have a nice day :-) from rank "<<my_rank<<endl;
    cout<<"--------------------------------------------"<<endl<<endl;
    */









// PART THEO
    // we need a proper definition for these things
    int NX = precs;
    int NY = UPP[my_rank] / precs;


    //variable declaration
    const int fullSize = NX * NY;
    std::vector<std::vector <double>> solutionU(2, std::vector<double>(fullSize,0));
    std::vector<double> ghostValues(fullSize,0);
    std::vector<double> ghostInNorth(NX,0), ghostInSouth(NX,0); // are dimensional allocations correct here?
    std::vector<double> ghostOutNorth(NX,0), ghostOutSouth(NX,0);
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
        MPI_Irecv( &ghostInNorth, NY, MPI_DOUBLE, idNorth, counter, comm1D, &requestNorth);
        MPI_Irecv( &ghostInSouth, NY, MPI_DOUBLE, idSouth, counter, comm1D, &requestSouth);

        // initiate send
        MPI_Send(&ghostOutNorth, NX,MPI_DOUBLE, idNorth, counter, comm1D);
        MPI_Send(&ghostOutSouth, NX,MPI_DOUBLE, idSouth, counter, comm1D);

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
                sum += A[i+NX*j] * solutionU[(counter+1)%2][i];
            }
            solutionU[counter%2][i] = (b[i] + ghostValues[i] - sum)/A[i+NX*i];
        }

        //calc runtime
        procRuntime += std::chrono::steady_clock::now() - start;
    }

    //if(counter == iterations) std::cout << "problem with counter not being interations" << std::endl;

    // prepare output
    std::vector<double> finalSolution(fullSize);
    std::vector<double> rhs(fullSize);
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

    the way NX and NY are defined:
    int NX = precs;
    int NY = UPP[my_rank] / precs;



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
    int root = 0;
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

//this is a test commit