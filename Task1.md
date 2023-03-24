# Task 1

## 1. Describe the advantages/disadvantages of a two-dimensional decomposition (compared to a one-dimensional decomposition).

The scaling behavior between 1D decomposition 2D decomposition is significantly different, as the cutting surface of the 2D subdomains is smaller than in the 1D case. But this will not affect the speedup much for small numbers of CPU cores but for large numbers. 

Example: 

take 250 resolution, 100 timesteps: 
t_s = time of communication of stecil
t_c = time of communication to other cell 
--> lets assume both take 1 time unit

T_p = parallel runtime 
T_s = 250^2 * t_s * 100 = runtime--> speedup: S = T_s/T_p 

1D Decomposition:
- 5 MPI-processes --> T_p = (250*200*t_s + 250 * 2 * t_c) * 100 --> S = 4.80
- 25 MPI-processes --> S = 20.83
- 125 MPI-processes --> S = 62.5
- 250 MPI processes --> S = 83.33

2D Decomposition:
- 5 MPI-processes --> S = 4.97
- 25 MPI-processes --> S = 24.38
- 125 MPI-processes --> S = 110.96
- 250 MPI processes --> S = 199.53

25 Processes:
 _
| | 
| |                       ____
| |   --> 2x 250         |    |  --> 4x 25
| |       = 500          |____|      = 1000
| |
|_|

125 Processes:
 _
| | 
| |                       ____
| |  --> 2x 250          |    |  --> 4x 2 
| |      = 500           |____|      = 4
| |
|_|

https://alexander.vondrous.de/?p=7

## 2. Discuss if the decomposition of the domain changes the order of computations performed during a single Jacobi iteration (i.e., if you expect a numerically identical result after each iteration, or not).



## 3. A generalization of the ghost layer approach would be to set the width of the ghost layer that is exchanged as a parameter W of the decomposition. This allows to perform W independent iterations before a communication of the ghost layers has to happen. Comment in which situation (w.r.t the available bandwidth or latency between MPI-processes) multiple independent iterations are potentially advantageous.



## 4. Assume a ghost layer with width W=1 (this is what you will later implement) and discuss if a data exchange between parts of the domain which are "diagonal neighbors" is required assuming a "5-point star-shaped stencil".

nicht direkt glaube ich, wenn bei 2d indirekt die Werte am jeweiligen eck durch horizontal und dann vertical bekommt diagonal element einen wert kommuniziert. Dieser würden den Wert aber eigentlich nicht brauchen. 

![](task_description/images/unitsquare_decomposition_1D_2D.png) 

## 5. How big is the sum of all L2 caches for 2 nodes of the IUE-cluster 

L2 (Level 2) cache is slower than the L1 cache but bigger in size. Where an L1 cache may measure in kilobytes, modern L2 memory caches measure in megabytes.
https://www.makeuseof.com/tag/what-is-cpu-cache/

Regular compute node (10x): 
2x INTEL Xeon Gold 6248, 2.5GHz, 20C/40T
Fat compute node (2x): 
2x INTEL Xeon Gold 6248, 2.5GHz, 20C/40T
https://www.iue.tuwien.ac.at/research/computing-infrastructure/

--> ein INTEL Xeon Gold 6248, 2.5GHz hat L2 cache von 16.0 MB
https://www.cpubenchmark.net/cpu.php?cpu=Intel+Xeon+Gold+6248+%40+2.50GHz&id=3517

Also for 2 nodes und je 2 INTEL = 4*16 = 64MB