#include <iostream>
#include <mpi.h> // https://www.open-mpi.org/doc/current/

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    int num_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::cout << "printing to stdout from rank=" << rank << "."<< std::endl;
    std::cerr << "printing to stdout from rank=" << rank << "."<< std::endl;

    MPI_Finalize();

  return 0;
}
