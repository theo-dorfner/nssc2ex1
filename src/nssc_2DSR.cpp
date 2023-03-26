#include<mpi.h>
#include"../include/functions_2D.h"
#include<iostream>

int main(int argc, char* argv[])
{
    int my_rank;
    int proc;
    MPI_Comm cart_comm;

    //A. Cart Creation
    //A.1 Initialization of MPI
    MPI_Init(&argc, &argv);

    string dimension = argv[1];
    int resolution = atoi(argv[2]);

    int precs = resolution - 2;
    //int iterations = atoi(argv[3]);

    MPI_Comm_size(MPI_COMM_WORLD, &proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    //A.2 Cart Creation
    vector<int> DIV;
    DIV = BiggestDivisors(proc);
    
    if(dimension == "1D" || (DIV[1] == 1 && DIV[0] == proc)) 
    //1D case or 2D case
    {
        //1D case
        std::cout<<"Ja moin wir sind im 1D case"<<std::endl;
        vectorprinter(DIV);
    } 
    else
    {
        //2D case
        int dim[2] = {DIV[0], DIV[1]};
        int periodical[2] = {1, 1};
        int reorder = 0;
        int coord[2];

        MPI_Cart_create(MPI_COMM_WORLD, 2, dim, periodical, reorder, &cart_comm);
        MPI_Cart_coords(cart_comm, my_rank, 2, coord);
    
        //3. Calculation of support functions

        vector<int>NX;
        //every dimension of each block in x direction
        vector<int>NY;
        //every dimension of each block in y direction
        vector<int>x_begin;
        //gives back index, where block begins in x_direction
        vector<int>y_begin;
        //gives back index, where block begins in y_direction
        vector<int>UPP;
        //Unknowns per Process (Block)
        double h;

        NX = NX_fun(DIV[0], precs); 
        NY = NY_fun(DIV[1], precs);
        x_begin = XBegin(NX, precs);
        y_begin = YBegin(NY, precs);
        UPP = UPP_fun(NX, NY);
        h = H(resolution);

        //4. b and A and u initialization
        int N = UPP[my_rank];
        vector<double>A(N*N, 0);
        vector<double>u(N, 0);
        vector<double>b;

        //std::cout << h << std::endl;
        A = Initialize_A0(A, N, NX[coord[1]], h);
        b = Initialize_b0(b, coord[1], coord[0], x_begin, y_begin, h);

        //print function for every rank
        if(my_rank==3){
        std::cout << N << "   " << NX[coord[1]] << std::endl;
        //std::cout << "my_rank = " << my_rank << "; my coords: " << coord[0] << coord[1] << std::endl;
        matrix_printer(A, N);
        //vectorprinter_d(b);
        }     
    }
    
    MPI_Finalize();
    
    //b check


    return 0;
}

