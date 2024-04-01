# Parallel Genetic Algorithm for the Traveling Salesman Problem (TSP)

## Overview
This project presents an innovative solution to the Traveling Salesman Problem (TSP) using a genetic algorithm enhanced with parallel computing capabilities via the Message Passing Interface (MPI). By integrating AI techniques with parallel processing, this approach significantly reduces computation time while exploring a diverse set of potential solutions.

## Features
- **Genetic Algorithm**: Employs evolutionary strategies to efficiently navigate the search space and identify optimal or near-optimal solutions to the TSP.
- **Parallel Processing**: Utilizes MPI to distribute computational tasks across multiple processors, enabling faster execution and scalability.
- **Dynamic Load Balancing**: Ensures efficient utilization of computational resources by balancing the workload among processors.
- **Custom MPI Functions**: Implements advanced MPI operations to enhance data communication and synchronization among processes.

## How It Works
1. **Initialization**: The program starts by reading the coordinates of cities from an input file and initializing a population of random routes.
2. **Fitness Evaluation**: The fitness of each route (or solution) is evaluated based on the total distance traveled.
3. **Genetic Operations**: Selection, crossover, and mutation operations are applied to generate new offspring, simulating the evolutionary process.
4. **Parallel Execution**: The population is distributed across multiple processes, where each process executes the genetic algorithm independently.
5. **Global Optimization**: Processes communicate via MPI to share their best-found solutions, collectively moving towards an optimal or near-optimal solution.

## Installation and Usage
1. Install MPI on your system (e.g., [MPICH](https://www.mpich.org/) or [OpenMPI](https://www.open-mpi.org/)).
2. Clone this repository:
3. Compile the source code using an MPI compiler (e.g., `mpicc`):  mpicxx -o parallel_tsp parallel_tsp.cpp
4. Run the program with MPI: mpirun -np <number_of_processes> ./parallel_tsp

## Innovation and Creativity
- **AI-Driven Approach**: The genetic algorithm, an AI technique, dynamically evolves solutions, efficiently exploring the vast search space of the TSP.
- **Parallelism with MPI**: Parallel execution with MPI not only accelerates the computation but also facilitates a broader exploration of solutions, increasing the chances of finding the global optimum.
- **Custom MPI Operations**: Tailored MPI functions enhance the algorithm's performance and accuracy by optimizing data exchange and synchronization processes.

## Contributing
We welcome contributions! If you have suggestions for improvements or want to contribute code, please feel free to submit a pull request or open an issue.

## License
This project is licensed under the [MIT License](https://opensource.org/license/MIT)- see the LICENSE file for details.
