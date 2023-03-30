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


int main(int argc, char* argv[])
{
    int my_rank;
    int proc;
    MPI_Comm comm1D;

    //A. Cart Creation
    //A.1 Initialization of MPI
    MPI_Init(&argc, &argv);

    string dimension = argv[1];
    int resolution = atoi(argv[2]);
    int iterations = atoi(argv[3]);

    int precs = resolution - 2;
    //int iterations = atoi(argv[3]);

    MPI_Comm_size(MPI_COMM_WORLD, &proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


    //A.2 Cart Creation
    
    if(dimension == "2D")
    {
        vector<int> DIV;
        DIV = BiggestDivisors(proc);
        //vectorprinter(DIV);

        if(DIV[0] == proc && DIV[1] == 1)
        {
            dimension = "1D";
        }
    }

//------------------------------------------1D CASE---------------------------------------//
    
    if(dimension == "1D") 
    //1D case or 2D case
    {
        //1D case
        if(my_rank == 1) std::cout<<"This is the 1D case!"<<std::endl;
        
        int ndims = 1;
        int dims[1];
        int wrap_around[1];
        int reorder;
        int ierr;

        //A. MPI_Cart creation

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
        //Check if resolution high enough for PROCS
        if(precs < proc)
        {
            if(my_rank == 0) cout << "Please choose a resolution that for N procs is at least N+2" << endl;
            return 0;
        }

        //B.1. Definition of Matrix A, vector u and vector b; size of arrays corresponds to UPP of the rank;
        std::vector<double>b;
        std::vector<double>u(UPP[my_rank],0);
        //starting value for vector u is all zero

        //B.2. Helpfunctions for the Initialization of A and b;
        double h = H(resolution);
        //calculating h
        vector<int>y_begin;
        //vector for the initialization of b
        y_begin = UPPtoYBegin(UPP, precs, proc);

        //B.3. Initialization of A and b;
        //A = Initialize_A0(A, UPP[my_rank], precs, h);
        b = Initialize_b0_1D(b, y_begin, precs, h, my_rank, proc);
        double alpha = 4 + 4 * M_PI * M_PI * h * h;

        // PART THEO
        // definition of matrix sice if unknown variables
        int NX = precs;
        int NY = UPP[my_rank] / precs;

        //variable declaration
        const int fullSize = NX * NY;
        std::vector<std::vector <double>> solutionU(2, std::vector<double>(fullSize,0));
        std::vector<double> ghostValues(fullSize,0);
        std::vector<double> ghostInNorth(NX,0), ghostInSouth(NX,0); // are dimensional allocations correct here?
        std::vector<double> ghostOutNorth(NX,0), ghostOutSouth(NX,0);
        int idNorth, idSouth, procID{my_rank};
        std::chrono::duration<double,std::nano> procRuntime{0};
        MPI_Request requestNorth;
        MPI_Request requestSouth;
        MPI_Status statusNorth;
        MPI_Status statusSouth;
        //int* collector[2] = {&idNorth, &idSouth};
        int startSouth = (NY-1)*NX;

        //MPI_Cart_coords(comm1D, my_rank, 2, &coords);
        //std::cout << procID << " " << coords << std::endl;
    
        // collect neighbour IDs
        //int MPI_Cart_shift(MPI_Comm comm_cart, int direction, int disp, int *rank_source, int *rank_dest)
        MPI_Cart_shift(comm1D, 1, -1, &procID, &idNorth); //vertical - north
        MPI_Cart_shift(comm1D, 1, +1, &procID, &idSouth); //vertical - south

        //for(auto &elem : collector)if(elem < 0)elem = &MPI_PROC_NULL;
        //std::cout << "on " << my_rank << " going north is " << idNorth << std::endl;
        //std::cout << "on " << my_rank << " going south is " << idSouth << std::endl;
        if(my_rank == 0){
            printf("jacobiMPI | resolution: %i; iterations: %i; dimension: %i; processes: %i",resolution,iterations,ndims,proc);
            std::cout << std::endl;
        }

        // start iterations
        for(int counter = 0; counter < iterations; ++counter){
            auto start = std::chrono::steady_clock::now(); // start runtime timing

            // generate outgoing ghost layers
        for(int i=0; i<NX;++i) ghostOutSouth[i] = solutionU[(counter+1)%2][startSouth + i];
        for(int i=0; i<NX;++i) ghostOutNorth[i] = solutionU[(counter+1)%2][i];

        // initiate non-blocking receive
        //MPI_Irecv( buf, count, datatype, source, tag, comm, [OUT] &request_handle);
        MPI_Irecv( &ghostInNorth[0], NX, MPI_DOUBLE, idNorth, counter, comm1D, &requestNorth);
        MPI_Irecv( &ghostInSouth[0], NX, MPI_DOUBLE, idSouth, counter, comm1D, &requestSouth);

        // initiate send
        MPI_Send(&ghostOutNorth[0], NX,MPI_DOUBLE, idNorth, counter, comm1D);
        MPI_Send(&ghostOutSouth[0], NX,MPI_DOUBLE, idSouth, counter, comm1D);

        // wait for receive
        MPI_Wait(&requestNorth,&statusNorth);
        MPI_Wait(&requestSouth,&statusSouth);

        // --- start jacobi-calculation
        // resolve neighbouring arrays into fullSize index position
        for(int i=0; i<NX;++i) ghostValues[startSouth + i] = (+1) * ghostInSouth[i]; //need to clarify +1/-1 here
        for(int i=0; i<NX;++i) ghostValues[i] = (+1) * ghostInNorth[i];

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

            solutionU[counter%2][i] = (b[i] + ghostValues[i]*h*h - sum)/alpha;
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
            //if(my_rank == 0) std::cout << finalSolution[i] << std::endl;
            rhs[i] = b[i] + ghostValues[i]*h*h;
        }
        // calculate mean runtime in seconds
        double meanRuntime = procRuntime.count()/(iterations*1.0*std::pow(10,9));


        // PART STEFFI

        std::vector <double> solution;                // initialise correct solution
        solution = Initialize_up(solution, y_begin, precs, h, my_rank, proc);

        // Initializations
        std::vector <double> residual_elements(NX*NY, 0);
        std::vector <double> error_elemets(NX*NY, 0);

        // Calculate Residual and Error
        residual_elements = Residual_Calc(NX, NY, resolution, finalSolution, rhs);
        error_elemets = Error_Calc(UPP[my_rank], finalSolution, solution);

        // Calculating Error and Residual Norm per Process
        double residualNorm_proc = NormL2_2(residual_elements);
        double residualMax_proc = NormInf(residual_elements);
        double errorNorm_proc = NormL2_2(error_elemets);
        double errorMax_proc = NormInf(error_elemets);

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
            double mean_runtime = runtime_sum/proc;
            auto residualNorm = sqrt(residualNorm_2);
            auto errorNorm = sqrt(errorNorm_2);

            // Output the result
            std::cout << std::scientific << "|residual|= " << residualNorm << std::endl;
            std::cout << std::scientific << "|residualMax|= " << residualMax << std::endl;
            std::cout << std::scientific << "|error|= " << errorNorm << std::endl;
            std::cout << std::scientific << "|errorMax|= " << errorMax << std::endl;
            std::cout << std::scientific << "average_runtime_per_iteration= " << mean_runtime << std::endl;
            std::cout << "---" << std::endl;

        }

    } 

//------------------------------------------2D CASE---------------------------------------//

    else if(dimension == "2D")
    {
        //2D case 
        if(my_rank == 0) std::cout<<"This is the 2D case"<<std::endl;

        vector<int> DIV;
        DIV = BiggestDivisors(proc);
    
        //2D case
        int dim[2] = {DIV[0], DIV[1]};
        int periodical[2] = {0, 0};
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
        double h;

        NX_rank = N_DIST(DIV[0], precs); 
        NY_rank = N_DIST(DIV[1], precs);
        x_begin = XBegin(NX_rank, precs);
        y_begin = YBegin(NY_rank, precs);
        int UPP = NX_rank[coord[1]] * NY_rank[coord[0]];
        h = H(resolution); 

        /*
        if(my_rank==1)
        {
            vectorprinter(NX_rank);
            vectorprinter(NY_rank);
            std::cout<<coord[0]<<" "<<coord[1]<<endl;
        }*/

        //std::cout<<my_rank<<" "<<NX_rank[coord[1]]<<NY_rank[coord[0]]<<" "<<UPP<<std::endl;
        vector<double>b;

        double alpha = 4 + 4 * M_PI * M_PI * h * h; 
        b = Initialize_b0_2D(b, coord[1], coord[0], x_begin, y_begin, h);

        // PART THEO
        // definition of matrix sice if unknown variables
        int NX = NX_rank[coord[1]];
        int NY = NY_rank[coord[0]];
        nprintf("my rank: %i; my coords: %i %i; nx, ny: %i %i ; my unknowns: %i \n",my_rank,coord[0],coord[1],NX,NY,UPP);

        
        //variable declaration
        std::vector<std::vector <double>> solutionU(2, std::vector<double>(UPP,0));
        std::vector<double> ghostValues(UPP,0);
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
        //printf("%i started carting\n",my_rank);
        //int MPI_Cart_shift(MPI_Comm comm_cart, int direction, int disp, int *rank_source, int *rank_dest)
        MPI_Cart_shift(comm1D, 1, -1, &procID, &idNorth); //vertical - north
        //printf("%i finished north %i\n",my_rank,idNorth);
        MPI_Cart_shift(comm1D, 1, +1, &procID, &idSouth); //vertical - south
        //printf("%i finished south %i\n",my_rank,idSouth);
        MPI_Cart_shift(comm1D, 0, -1, &procID, &idWest); //horizontal - west
        //printf("%i finished west %i\n",my_rank,idWest);
        MPI_Cart_shift(comm1D, 0, +1, &procID, &idEast); //horizontal - east
        //printf("%i finished east %i\n",my_rank,idEast);

        //for(auto &elem : collector)if(elem < 0)elem = &MPI_PROC_NULL;

        //std::cout << "on " << my_rank << " going north is " << idNorth << std::endl;
        //std::cout << "on " << my_rank << " going south is " << idSouth << std::endl;
        if(my_rank == 0)std::cout << printf("jacobiMPI | resolution: %i; iterations: %i; processes: %i",resolution,iterations,proc) << std::endl;

        
        // start iterations
        for(int counter = 0; counter < iterations; ++counter){
            auto start = std::chrono::steady_clock::now(); // start runtime timing
        
        
        // generate outgoing ghost layers
        for(int i=0; i<NX;++i) ghostOutSouth[i] = solutionU[(counter+1)%2][startSouth + i];
        for(int i=0; i<NX;++i) ghostOutNorth[i] = solutionU[(counter+1)%2][i];
        for(int i=0; i<NY;++i) ghostOutWest[i] = solutionU[(counter+1)%2][i*NX];
        for(int i=1; i<NY+1;++i) ghostOutEast[i] = solutionU[(counter+1)%2][NX*i - 1];

        //printf("%i starts receiving\n",my_rank);

        // initiate non-blocking receive
        //MPI_Irecv( buf, count, datatype, source, tag, comm, [OUT] &request_handle);
        //printf("%i start receiving north\n",my_rank);
        MPI_Irecv( &ghostInNorth[0], NX, MPI_DOUBLE, idNorth, counter, comm1D, &requestNorth); // still needs changing of comm1D
        //printf("%i start receiving south\n",my_rank);
        MPI_Irecv( &ghostInSouth[0], NX, MPI_DOUBLE, idSouth, counter, comm1D, &requestSouth);
        //printf("%i start receiving west\n",my_rank);
        MPI_Irecv( &ghostInWest[0], NY, MPI_DOUBLE, idWest, counter, comm1D, &requestWest);
        //printf("%i start receiving east\n",my_rank);
        MPI_Irecv( &ghostInEast[0], NY, MPI_DOUBLE, idEast, counter, comm1D, &requestEast);

        //printf("%i passed receives\n",my_rank);

        
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
        

        //printf("%i passed waits\n",my_rank);

        //printf("proc %i; cnt %i; status north: ",my_rank,counter);
        //std::cout << statusNorth.MPI_ERROR << std::endl;

        // --- start jacobi-calculation
        // resolve neighbouring arrays into fullSize index position
        for(int i=0; i<NX;++i) ghostValues[startSouth + i] = (+1) * ghostInSouth[i]; //need to clarify +1/-1 here
        for(int i=0; i<NX;++i) ghostValues[i] = (+1) * ghostInNorth[i];
        for(int i=0; i<NY;++i) ghostValues[i*NX] = (+1) * ghostInWest[i];
        for(int i=1; i<NY+1;++i) ghostValues[NX*i - 1] = (+1) * ghostInEast[i];

        int precs = NX;

        

        // actually calculate jacobi
        // das sind andere i,j als in den grid-koordinaten (diese hier sind nur intern für berechnungen der Matrix)
        for(int i=0; i<UPP; i++){
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
            else if(i>=precs && i<(UPP-precs))
            {
                sum +=  - solutionU[(counter+1)%2][i-precs]
                        - solutionU[(counter+1)%2][i+precs];
                
                if((i+1)%precs != 0)    sum += - solutionU[(counter+1)%2][i+1];
                if(i%precs != 0)        sum += - solutionU[(counter+1)%2][i-1];
            }
            else if(i>=(UPP-precs) && i<UPP-1)
            {
                sum +=  - solutionU[(counter+1)%2][i-precs];
                
                if((i+1)%precs != 0)    sum += - solutionU[(counter+1)%2][i+1];
                if(i%precs != 0)        sum += - solutionU[(counter+1)%2][i-1];
            }
            else if(i==UPP-1)
            {
                sum +=  - solutionU[(counter+1)%2][UPP-2] 
                        - solutionU[(counter+1)%2][UPP-1-precs];
            }
            
            //std::cout << ghostValues[i] << std::endl;
            solutionU[counter%2][i] = (b[i] + ghostValues[i]*h*h - sum)/alpha;
            }

            //calc runtime
            procRuntime += std::chrono::steady_clock::now() - start;
          
        } 
        
    
        //if(counter == iterations) std::cout << "problem with counter not being interations" << std::endl;

        // prepare output
        std::vector<double> finalSolution(UPP,0);
        std::vector<double> rhs(UPP,0);
        for(int i=0; i < UPP; ++i){
            finalSolution[i] = solutionU[iterations % 2][i];
            //if(my_rank == 0) std::cout << finalSolution[i] << std::endl;
            rhs[i] = b[i] + ghostValues[i]*h*h;
        }
        // calculate mean runtime in seconds
        double meanRuntime = procRuntime.count()/(std::pow(10,9)*iterations*1.0);

        

        // PART STEFFI

        

        std::vector <double> solution;            // initialise correct solution
        solution = Initialize_up(solution, y_begin, precs, h, my_rank, proc);

        

        // Initializations
        std::vector <double> residual_elements(NX*NY, 0);
        std::vector <double> error_elemets;

        
        // Calculate Residual and Error
        residual_elements = Residual_Calc(NX, NY, resolution, finalSolution, rhs);
        //error_elemets = Error_Calc(UPP, finalSolution, solution); ATTENTION!

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
        double mean_runtime = runtime_sum/proc;
        auto residualNorm = sqrt(residualNorm_2);
        //auto errorNorm = sqrt(errorNorm_2);

        // Output the result
        std::cout << std::endl;
        std::cout << std::scientific << "|residual|= " << residualNorm << std::endl;
        std::cout << std::scientific << "|residualMax|= " << residualMax << std::endl;
        //std::cout << std::scientific << "|error|= " << errorNorm << std::endl;
        //std::cout << std::scientific << "|errorMax|= " << errorMax << std::endl;
        std::cout << std::scientific << "average_runtime_per_iteration= " << mean_runtime << std::endl;

        
        }
    }
    

    MPI_Finalize();

}