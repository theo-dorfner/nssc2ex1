#include <stdio.h>
#include <mpi.h>
#include <iostream>

#include <cmath>
#include <vector>
#include <array>
#include <limits>
#include <chrono>

#include "../include/functions_2D.h"
#include "../include/functions_sn.h"
#include "../include/jacobi.hpp"

using matrixType = std::vector <std::vector <double>>;

using Clock = std::chrono::system_clock;

int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    int my_rank;
    int procs;
    MPI_Comm comm1D;

    //A. Cart Creation
    //A.1 Initialization of MPI

    string dimension = argv[1];
    int resolution = atoi(argv[2]);
    int iterations = atoi(argv[3]);

    int precs = resolution - 2;
    //int iterations = atoi(argv[3]);

    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    //A.2 Cart Creation
    vector<int> DIV;
    DIV = BiggestDivisors(procs);
    
    //2D case
    int dim[2] = {DIV[0], DIV[1]};
    int periodical[2] = {1, 1};
    int reorder = 0;
    int coord[2];

    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, periodical, reorder, &comm1D);
    MPI_Cart_coords(comm1D, my_rank, 2, coord);
    
    //3. Calculation of support functions

    vector<int>NX_rank;
    //every dimension of each block in x direction
    vector<int>NY_rank;
    //every dimension of each block in y direction
    vector<int>x_begin;
    //gives back index, where block begins in x_direction
    vector<int>y_begin;
    //gives back index, where block begins in y_direction
    vector<int>UPP;
    //Unknowns per Process (Block)
    double h;

    NX_rank = NX_fun(DIV[0], precs); 
    NY_rank = NY_fun(DIV[1], precs);
    x_begin = XBegin(NX_rank, precs);
    y_begin = YBegin(NY_rank, precs);
    UPP = UPP_fun(NX_rank, NY_rank);
    h = H(resolution);

    //4. b and A and u initialization
    //int N = UPP[my_rank];
    //vector<double>A(N*N, 0);
    vector<double>u(UPP[my_rank], 0);
    vector<double>b;

    //std::cout << h << std::endl;
    //A = Initialize_A0(A, N, NX_rank[coord[1]], h);
    double alpha = 4 + 4 * M_PI * M_PI * h * h; 
    b = Initialize_b0(b, coord[1], coord[0], x_begin, y_begin, h);


// PART THEO
    // definition of matrix sice if unknown variables
    const int NX = NX_rank[my_rank];
    const int NY = NX_rank[my_rank];


    //variable declaration
    const int fullSize = NX * NY;
    std::vector<std::vector <double>> solutionU(2, std::vector<double>(fullSize,0));
    std::vector<double> ghostValues(fullSize,0);
    std::vector<double> ghostInNorth(NX,0), ghostInSouth(NX,0), ghostInWest(NY,0), ghostInEast(NY,0); // are dimensional allocations correct here?
    std::vector<double> ghostOutNorth(NX,0), ghostOutSouth(NX,0), ghostOutWest(NY,0), ghostOutEast(NY,0);
    int idNorth, idSouth, idWest, idEast, procID{my_rank};
    std::chrono::duration<double> procRuntime{0};
    MPI_Request requestNorth;
    MPI_Request requestSouth;
    MPI_Request requestWest;
    MPI_Request requestEast;
    MPI_Status statusNorth;
    MPI_Status statusSouth;
    MPI_Status statusWest;
    MPI_Status statusEast;
    
    //int* collector[2] = {&idNorth, &idSouth};
    int startSouth = (NY-1)*NX;

    //MPI_Cart_coords(comm1D, my_rank, 2, &coords);
    //std::cout << procID << " " << coords << std::endl;
    
    // collect neighbour IDs
    //int MPI_Cart_shift(MPI_Comm comm_cart, int direction, int disp, int *rank_source, int *rank_dest)
    MPI_Cart_shift(comm1D, 1, -1, &procID, &idNorth); //vertical - north
    MPI_Cart_shift(comm1D, 1, +1, &procID, &idSouth); //vertical - south
    MPI_Cart_shift(comm1D, 0, -1, &procID, &idWest); //horizontal - west
    MPI_Cart_shift(comm1D, 0, +1, &procID, &idEast); //horizontal - east

    //for(auto &elem : collector)if(elem < 0)elem = &MPI_PROC_NULL;

    //std::cout << "on " << my_rank << " going north is " << idNorth << std::endl;
    //std::cout << "on " << my_rank << " going south is " << idSouth << std::endl;
    if(my_rank == 0)std::cout << printf("jacobiMPI | resolution: %i; iterations: %i; processes: %i",resolution,iterations,procs) << std::endl;

    
    // start iterations
    for(int counter = 0; counter < iterations; ++counter){
        auto start = std::chrono::steady_clock::now(); // start runtime timing

        // generate outgoing ghost layers
        for(int i=0; i<NX;++i) ghostOutSouth[i] = solutionU[(counter+1)%2][startSouth + i];
        for(int i=0; i<NX;++i) ghostOutNorth[i] = solutionU[(counter+1)%2][i];
        for(int i=0; i<NY;++i) ghostOutWest[i] = solutionU[(counter+1)%2][i*NX];
        for(int i=1; i<NY+1;++i) ghostOutEast[i] = solutionU[(counter+1)%2][NX*i - 1];

        // initiate non-blocking receive
        //MPI_Irecv( buf, count, datatype, source, tag, comm, [OUT] &request_handle);
        MPI_Irecv( &ghostInNorth[0], NX, MPI_DOUBLE, idNorth, counter, comm1D, &requestNorth); // still needs changing of comm1D
        MPI_Irecv( &ghostInSouth[0], NX, MPI_DOUBLE, idSouth, counter, comm1D, &requestSouth);
        MPI_Irecv( &ghostInWest[0], NY, MPI_DOUBLE, idWest, counter, comm1D, &requestWest);
        MPI_Irecv( &ghostInEast[0], NY, MPI_DOUBLE, idEast, counter, comm1D, &requestEast);

        // initiate send
        MPI_Send(&ghostOutNorth[0], NX,MPI_DOUBLE, idNorth, counter, comm1D);
        MPI_Send(&ghostOutSouth[0], NX,MPI_DOUBLE, idSouth, counter, comm1D);
        MPI_Send(&ghostOutWest[0], NY,MPI_DOUBLE, idWest, counter, comm1D);
        MPI_Send(&ghostOutEast[0], NY,MPI_DOUBLE, idEast, counter, comm1D);

        
        // wait for receive
        MPI_Wait(&requestNorth,&statusNorth);
        MPI_Wait(&requestSouth,&statusSouth);
        MPI_Wait(&requestWest,&statusWest);
        MPI_Wait(&requestEast,&statusEast);
        
        

        // --- start jacobi-calculation
        // resolve neighbouring arrays into fullSize index position
        for(int i=0; i<NX;++i) ghostValues[startSouth + i] = (+1) * ghostInSouth[i]; //need to clarify +1/-1 here
        for(int i=0; i<NX;++i) ghostValues[i] = (+1) * ghostInNorth[i];
        for(int i=0; i<NY;++i) ghostValues[i*NX] = (+1) * ghostInWest[i];
        for(int i=1; i<NY+1;++i) ghostValues[NX*i - 1] = (+1) * ghostInEast[i];

        int precs = NX;

        // actually calculate jacobi
        // das sind andere i,j als in den grid-koordinaten (diese hier sind nur intern für berechnungen der Matrix)
        for(int i=0; i<fullSize; i++){
            double sum{0};
            if(i==0)
            {
                sum +=  - solutionU[(counter+1)%2][1] 
                        - solutionU[(counter+1)%2][precs];
            }
            else if(i>0 && i<precs)
            {
                sum +=  - solutionU[(counter+1)%2][i+precs];           
                
                if((i+1)%precs != 0)    sum += - solutionU[(counter+1)%2][i+1];
                if(i%precs != 0)        sum += - solutionU[(counter+1)%2][i-1];
            }
            else if(i>=precs && i<(UPP[my_rank]-precs))
            {
                sum +=  - solutionU[(counter+1)%2][i-precs]
                        - solutionU[(counter+1)%2][i+precs];
                
                if((i+1)%precs != 0)    sum += - solutionU[(counter+1)%2][i+1];
                if(i%precs != 0)        sum += - solutionU[(counter+1)%2][i-1];
            }
            else if(i>=(UPP[my_rank]-precs) && i<UPP[my_rank]-1)
            {
                sum +=  - solutionU[(counter+1)%2][i-precs];
                
                if((i+1)%precs != 0)    sum += - solutionU[(counter+1)%2][i+1];
                if(i%precs != 0)        sum += - solutionU[(counter+1)%2][i-1];
            }
            else if(i==UPP[my_rank]-1)
            {
                sum +=  - solutionU[(counter+1)%2][UPP[my_rank]-2] 
                        - solutionU[(counter+1)%2][UPP[my_rank]-1-precs];
            }

            /*
            for(int j=0; j<fullSize; j++) {
                if(i==j) continue;
                sum += A_at(A,X,Y,i,j) * solutionU[(counter+1)%2][j];
            }
            */
            //std::cout << ghostValues[i] << std::endl;
            solutionU[counter%2][i] = (b[i] + ghostValues[i]*h*h - sum)/alpha;
        }

        /*
        // actually calculate jacobi
        // das sind andere i,j als in den grid-koordinaten (diese hier sind nur intern für berechnungen der Matrix)
        for(int i=0; i<fullSize; i++){
            double sum{0};
            for(int j=0; j<fullSize; j++) {
                if(i==j) continue;
                sum += A[j+fullSize*i] * solutionU[(counter+1)%2][j];
            }
            //std::cout << ghostValues[i] << std::endl;
            solutionU[counter%2][i] = (b[i] + ghostValues[i]*h*h - sum)/A[i+fullSize*i];
        }
        */

        //calc runtime
        procRuntime += std::chrono::steady_clock::now() - start;
    }

    //if(counter == iterations) std::cout << "problem with counter not being interations" << std::endl;



    // prepare output
    std::vector<double> finalSolution(fullSize);
    std::vector<double> rhs(fullSize);
    for(int i=0; i < fullSize; ++i){
        finalSolution[i] = solutionU[iterations % 2][i];
        //if(my_rank == 0) std::cout << finalSolution[i] << std::endl;
        rhs[i] = b[i] + ghostValues[i]*h*h;
    }
    // calculate mean runtime in seconds
    double meanRuntime = procRuntime.count()/(std::pow(10,9)*iterations*1.0);


    


// PART STEFFI

    std::vector <double> solution;            // initialise correct solution
    solution = Initialize_up(solution, y_begin, precs, h, my_rank, procs);

    // Initializations
    std::vector <double> residual_elements(NX*NY, 0);
    std::vector <double> error_elemets(NX*NY, 0);

    // Calculate Residual and Error
    residual_elements = Residual_Calc(NX, NY, resolution, finalSolution, rhs);
    error_elemets = Error_Calc(NX, NY, finalSolution, solution);

    // Clculating Error and Residual Norm per Process
    double residualNorm_proc = NormL2_2(residual_elements);
    double residualMax_proc = NormInf(residual_elements);
    double errorNorm_proc = NormL2_2(error_elemets);
    double errorMax_proc = NormInf(error_elemets);

    //std::cout << my_rank << " : "<< errorMax_proc << std::endl;
    //if(my_rank == 6) matrix_printer(A, UPP[my_rank]);

    // Initialize for combining
    int root = 0;
    double residualNorm_2;
    double residualMax;
    double errorNorm_2;
    double errorMax;
    double runtime_sum;

   //if(my_rank == 0)vector_printer(error_elemets);

    // gather all error and residual and runtime to combine and sum up 
    MPI_Reduce( &residualNorm_proc, &residualNorm_2, 1, MPI_DOUBLE,
                MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce( &residualMax_proc, &residualMax, 1, MPI_DOUBLE,
                MPI_MAX, 0, MPI_COMM_WORLD);

    MPI_Reduce( &errorNorm_proc, &errorNorm_2, 1, MPI_DOUBLE,
                MPI_SUM, root, MPI_COMM_WORLD);
    
    MPI_Reduce( &errorMax_proc, &errorMax, 1, MPI_DOUBLE,
                MPI_MAX, 0, MPI_COMM_WORLD);

    MPI_Reduce( &meanRuntime, &runtime_sum, 1, MPI_DOUBLE,
                MPI_SUM, root, MPI_COMM_WORLD);
        
    if(my_rank == 0){

        // Calculating NormL2 and Infinite Norm of Residual and Error
        double mean_runtime = runtime_sum/procs;
        auto residualNorm = sqrt(residualNorm_2);
        auto errorNorm = sqrt(errorNorm_2);

        // Output the result
        std::cout << std::scientific << "|residual|= " << residualNorm << std::endl;
        std::cout << std::scientific << "|residualMax|= " << residualMax << std::endl;
        std::cout << std::scientific << "|error|= " << errorNorm << std::endl;
        std::cout << std::scientific << "|errorMax|= " << errorMax << std::endl;
        std::cout << std::scientific << "average_runtime_per_iteration= " << mean_runtime << std::endl;
    
    }
    

    MPI_Finalize();

    

}
