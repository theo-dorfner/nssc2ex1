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

    //A. MPI_Cart creation

    //A.1. Reading in values;
    int resolution = atoi(argv[1]);
    int precs = resolution - 2; //This is N_x
    int Unknowns = precs * precs;
    //int iterations = atoi(argv[2]);


    //A.2. Reading out proc and rank; create cart and check if successful
    MPI_Comm_size(MPI_COMM_WORLD, &proc);
    //saves number of procs on "proc"
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    //saves rank of MPI on "my_rank"

    dims[0] = 0;

    MPI_Dims_create(proc, ndims, dims);

    wrap_around[0] = 0;
    reorder = 1;
    ierr = 0;
    ierr = MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, wrap_around, reorder, &comm1D);
    //creating MPI cart with one dimension

    if(ierr != 0) printf("ERROR[%i] creating CART\n", ierr);
    //returns error if cart creation wasn't successful


    // B: Definition of A and b

    vector<int>UPP;
    UPP = UnknownsPerProc(UPP, precs, proc);
    //calculates the number of Unknowns per Process (saves them in a vector
    int NY = UPP[my_rank] / precs;
    //Check if resolution high enough for procs; if not: cancel program!
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

    cout << "Hello, my rank is " << my_rank << " and the number of unknowns are " << UPP[my_rank] << endl << endl;
    cout << "This is my vector b: " << endl;
    vector_printer(b);
    cout <<"And this is my matrix A: " << endl;
    matrix_printer(A, UPP[my_rank]);
    cout<<"Have a nice day :-) from rank "<<my_rank<<endl;
    cout<<"--------------------------------------------"<<endl<<endl;

    MPI_Finalize();

}
