stenciljacobi: cppMain.o
	g++ -Wall -std=c++2a -pedantic -O2 -C -o stenciljacobi cppMain.o 
	make clean

cppMain.o: cppJacobi.cpp cppJacobi.hpp
	g++ -Wall -std=c++2a -pedantic -O2 -c -o cppMain.o cppJacobi.cpp 

clean:
	rm -f cppMain.o
