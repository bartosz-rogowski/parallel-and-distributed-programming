# Parallel and distributed systems

Projects for university course "Parallel and distributed systems".

## Table of contents
1. [Parallel MPI Application](#parallel)
1. [Distributed Application](#distributed)

---

## 1. Parallel MPI Application: Minimal Spanning Tree - Prim's Algorithm <a name="parallel"></a>

Prim's algorithm (Parallel Version) written within OpenMPI. 

How to run a program (via *Makefile*):
* `make prepare` - for preparing environment (adjusted to the laboratory conditions)
* `make` - for compiling the program
* `make run` - to run the default version of program, takes 3 args:
  - `n` - number of processes
  - `in` - path to the file with input data
  - `out` *optional* - path to the file with output data (adjacency matrix of MST)
* `make run0_` (`_` means number between 1 and 5)

[*link to documentation*](./parallel/doc.pdf)

---

## 2. Distributed <*type to be established*> Application: <*subject to be established*> <a name="distributed"></a>

Prim's algorithm (Distributed Version) written within UPC++. 

How to run a program (via *Makefile*):
* `make prepare` - for preparing environment (adjusted to the laboratory conditions)
* `make compile` - for compiling the program
* `make run` - to run the default version of program, takes 3 args:
  - `n` - number of processes
  - `in` - path to the file with input data
  - `out` *optional* - path to the file with output data (adjacency matrix of MST)
* `make run0_` (`_` means number between 1 and 5)


[*link to documentation*](./distributed/doc.pdf)
