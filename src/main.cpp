#include <assert.h>
#include <iomanip>
#include <iostream>
#include <vector>

#include "arguments.hpp"
#include "solver.hpp"
 
int main(int argc, char **argv) {

  // parse command line arguments
  auto resolution = convertTo<int>(1, 32, argc, argv);
  auto iterations = convertTo<int>(2, 800, argc, argv);

  std::cout << "resolution=" << resolution << std::endl;
  std::cout << "iterations=" << iterations << std::endl;

  assert(resolution > 0);
  assert(iterations > 0);

  PoissonJacobiStencil(resolution, iterations);

  return 0;
}
