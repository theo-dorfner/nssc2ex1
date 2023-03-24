#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include "../include/functions.h"
#include <vector>

using namespace std;

int main(int argc, char* argv[]) {

    int ndims = 1;
    int size;
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
    //int iterations = atoi(argv[2]);


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
    vector<int>UPP;
    UPP = UnknownsPerProc(UPP, precs, proc); //calculates how big dimensions of arrays are for each rank!

    //IMPLEMENT FUNCTION(PROOF IF ENOUGH PRECS FOR PROCS)

    vector<double>b;
    vector<double>u(UPP[my_rank],0);
    vector<double>A(UPP[my_rank]*UPP[my_rank],0); //this is a Matrix --> initialize with A[j*N+i]
    //double* A[UPP[my_rank]]; //matrix; initialize with A[j*N + i]

    //u = Initialize_u0(u, UPP[my_rank]); //for "normal" application
    //b = Initialize_b0(b, UPP[my_rank]);
    double h = H(resolution);

    A = Initialize_A0(A, UPP[my_rank], precs, h);

    vector<int>y_begin;
    y_begin = UPPtoYBegin(UPP, precs, proc);

    vector<double>Y_vals;
    b = Initialize_b0(b, y_begin, precs, h, my_rank, proc);

    //vector_printer_int(y_begin);

    cout<< "rank: " << my_rank << endl;
    vector_printer(b);
    cout << endl;


    //vector_printer(b);
    //matrix_printer(A, UPP[my_rank]);




    //b = Initialize_v_random(b, UPP[my_rank]); //for Steffis application
    //u = Initialize_v_random(u, UPP[my_rank]);
    //A = Initialize_A_random(A, UPP[my_rank]);

    /*
    cout << "Hello, my rank is " << my_rank << " and the number of unknowns are " << UPP[my_rank] << endl << endl;
    cout << "This is my vector b: " << endl;
    vector_printer(b);
    cout << "This is my vector u: " << endl;
    vector_printer(u);
    cout <<"And this is my matrix A: " << endl;
    matrix_printer(A, UPP[my_rank]);
    cout<<"Have a nice day :-) from rank "<<my_rank<<endl;
    cout<<"--------------------------------------------"<<endl<<endl;
    */


    MPI_Finalize();

}

/* "Hello world" for Cart; shows coordinates of each rank plus up and down "neighbour"
MPI_Cart_coords(comm1D, my_rank, ndims, &coord);

int up, down;
MPI_Cart_shift(comm1D, 0, 1, &up, &down);

std::cout<<"hello my rank is: "<< my_rank << " and my coordinate is: " << coord<<std::endl;
std::cout<<"rank up is: "<<up<<" rank down is: "<< down << std::endl << std::endl;
*/
