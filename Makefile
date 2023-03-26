jacobiMPI: ../src/jacobiMPI.cpp ../include/functions_sn.h ../include/functions.h ../include/jacobi.hpp
	mpic++ -std=c++17 -O3 -Wall -pedantic -march=native -ffast-math ../src/jacobiMPI.cpp -o jacobiMPI
