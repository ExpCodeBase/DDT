# DDT

This is the code base for distributed truss decomposition in directed graphs.

## File structure

```
/DDT
    /src
        main.cpp
        partitioner.cpp
        partitioner.h
        peeler.cpp
        peeler.h
        graph.cpp
        graph.h
        mpi_utils.cpp
        mpi_utils.h
    /data
        graph.e
    README.md
```

## Requirements

- C++11 or later
- A C++ compiler (e.g., g++, clang++)
- OpenMPI

## Setup

1. Clone the repository or download the source code.
2. Ensure the `graph.e` file is placed in the `data` directory.
3. Install OpenMPI and ensure it is configured correctly.

## Compilation

Navigate to the `src` directory and run the following command:

```bash
mpic++ -o main main.cpp graph.cpp partitioner.cpp peeler.cpp mpi_utils.cpp
```

## Running the Program
After compilation, run the executable using MPI:

```bash
mpirun -np <number_of_processes> ./main
```

Replace <number_of_processes> with the number of processes you want to run.

## Input File

The input file graph.e should be formatted such that each line contains two integers separated by a space, representing an edge from the first integer to the second integer.