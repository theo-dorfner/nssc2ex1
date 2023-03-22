# Task 1

## 1. Describe the advantages/disadvantages of a two-dimensional decomposition (compared to a one-dimensional decomposition).

more parallelization, symmetric matrices not rectangle maybe better?

## 2. Discuss if the decomposition of the domain changes the order of computations performed during a single Jacobi iteration (i.e., if you expect a numerically identical result after each iteration, or not).



## 3. A generalization of the ghost layer approach would be to set the width of the ghost layer that is exchanged as a parameter W of the decomposition. This allows to perform W independent iterations before a communication of the ghost layers has to happen. Comment in which situation (w.r.t the available bandwidth or latency between MPI-processes) multiple independent iterations are potentially advantageous.



## 4. Assume a ghost layer with width W=1 (this is what you will later implement) and discuss if a data exchange between parts of the domain which are "diagonal neighbors" is required assuming a "5-point star-shaped stencil".



## 5. How big is the sum of all L2 caches for 2 nodes of the IUE-cluster 

