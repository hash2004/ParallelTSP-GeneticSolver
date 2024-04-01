#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <random>
#include <cstring>
#include <mpi.h>

constexpr int MAX_CITIES = 100;
constexpr int MAX_POPULATION = 50;

struct City {
    double x, y;
};

double distances[MAX_CITIES][MAX_CITIES];
City cities[MAX_CITIES];
int population[MAX_POPULATION][MAX_CITIES];
int* populationPtrs[MAX_POPULATION];

double euclideanDistance(const City& a, const City& b) {
    return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

void calculateDistances(int cityCount) {
    for (int i = 0; i < cityCount; ++i) {
        for (int j = 0; j < cityCount; ++j) {
            distances[i][j] = euclideanDistance(cities[i], cities[j]);
        }
    }
}

double routeDistance(const int* route, int cityCount) {
    double totalDistance = 0.0;
    for (int i = 0; i < cityCount - 1; ++i) {
        totalDistance += distances[route[i]][route[i + 1]];
    }
    totalDistance += distances[route[cityCount - 1]][route[0]];  // Complete the loop
    return totalDistance;
}

void mutate(int* route, int cityCount) {
    std::random_device rd;
    std::mt19937 g(rd());
    std::uniform_int_distribution<> dis(0, cityCount - 1);
    int i = dis(g);
    int j = dis(g);
    std::swap(route[i], route[j]);
}

void crossover(const int* parent1, const int* parent2, int* child, int cityCount) {
    std::fill_n(child, cityCount, -1);

    std::random_device rd;
    std::mt19937 g(rd());
    std::uniform_int_distribution<> dis(0, cityCount - 1);

    int start = dis(g);
    int end = dis(g);
    if (start > end) std::swap(start, end);

    for (int i = start; i <= end; ++i) {
        child[i] = parent1[i];
    }

    int index = (end + 1) % cityCount;
    for (int i = 0; i < cityCount; ++i) {
        int city = parent2[i];
        if (std::find(child, child + cityCount, city) == child + cityCount) {
            while (child[index] != -1) {
                index = (index + 1) % cityCount;
            }
            child[index] = city;
        }
    }
}

void my_MPI_Allgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
                      void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // Each process sends its data to all other processes
    for (int i = 0; i < size; i++) {
        if (i == rank) {
            // Copy own data into the right position in recvbuf
            memcpy((char*)recvbuf + rank * recvcount * sizeof(sendtype), sendbuf, sendcount * sizeof(sendtype));
        } else {
            MPI_Send(sendbuf, sendcount, sendtype, i, 0, comm);
        }
    }

    // Receive data from all other processes
    for (int i = 0; i < size; i++) {
        if (i != rank) {
            MPI_Recv((char*)recvbuf + i * recvcount * sizeof(recvtype), recvcount, recvtype, i, 0, comm, MPI_STATUS_IGNORE);
        }
    }
}
void my_MPI_Allgatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
                       void* recvbuf, const int* recvcounts, const int* displs, 
                       MPI_Datatype recvtype, MPI_Comm comm) {
    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // Send own data to all other processes
    for (int i = 0; i < size; i++) {
        if (i == rank) {
            // Copy own data into the right position in recvbuf
            memcpy((char*)recvbuf + displs[rank] * sizeof(recvtype), sendbuf, sendcount * sizeof(sendtype));
        } else {
            MPI_Send(sendbuf, sendcount, sendtype, i, 0, comm);
        }
    }

    // Receive data from all other processes
    for (int i = 0; i < size; i++) {
        if (i != rank) {
            MPI_Recv((char*)recvbuf + displs[i] * sizeof(recvtype), recvcounts[i], recvtype, i, 0, comm, MPI_STATUS_IGNORE);
        }
    }
}
void my_MPI_Alltoall(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
                     void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    int sendtype_size, recvtype_size;
    MPI_Type_size(sendtype, &sendtype_size);
    MPI_Type_size(recvtype, &recvtype_size);

    // Send and receive data from all processes
    for (int i = 0; i < size; i++) {
        const char* sendPtr = (const char*)sendbuf + i * sendcount * sendtype_size;
        char* recvPtr = (char*)recvbuf + i * recvcount * recvtype_size;
        
        if (i == rank) {
            // Copy own data into the right position in recvbuf
            memcpy(recvPtr, sendPtr, sendcount * sendtype_size);
        } else {
            MPI_Send(sendPtr, sendcount, sendtype, i, 0, comm);
            MPI_Recv(recvPtr, recvcount, recvtype, i, 0, comm, MPI_STATUS_IGNORE);
        }
    }
}

void my_MPI_Alltoallv(const void* sendbuf, const int* sendcounts, const int* sdispls,
                      MPI_Datatype sendtype, void* recvbuf, const int* recvcounts,
                      const int* rdispls, MPI_Datatype recvtype, MPI_Comm comm) {
    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // Each process sends its data to every other process and receives from every other process
    for (int i = 0; i < size; i++) {
        char* sendPtr = (char*)sendbuf + sdispls[i] * sizeof(sendtype);
        MPI_Send(sendPtr, sendcounts[i], sendtype, i, 0, comm);

        char* recvPtr = (char*)recvbuf + rdispls[i] * sizeof(recvtype);
        MPI_Recv(recvPtr, recvcounts[i], recvtype, i, 0, comm, MPI_STATUS_IGNORE);
    }
}
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::ifstream file("input.txt");
    int cityCount = 0;

    if (rank == 0) {
        double x, y;
        while (file >> x >> y) {
            cities[cityCount] = {x, y};
            cityCount++;
        }
    }

    // Broadcast the number of cities to all processes
    MPI_Bcast(&cityCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Now distribute the city data to all processes
    MPI_Bcast(cities, cityCount * sizeof(City), MPI_BYTE, 0, MPI_COMM_WORLD);

    calculateDistances(cityCount);

    // Initialize populations
    for (int i = 0; i < MAX_POPULATION; ++i) {
        populationPtrs[i] = population[i];
        for (int j = 0; j < cityCount; ++j) {
            population[i][j] = j;
        }
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(population[i], population[i] + cityCount, g);
    }

    const int generations = 1000;
for (int i = 0; i < generations; ++i) {
    // Sort the population based on the route distance
    std::sort(populationPtrs, populationPtrs + MAX_POPULATION, [cityCount](const int* a, const int* b) {
        return routeDistance(a, cityCount) < routeDistance(b, cityCount);
    });

    int newPopulation[MAX_POPULATION][MAX_CITIES];
    for (int j = 0; j < MAX_POPULATION / 2; ++j) {
        int child1[MAX_CITIES];
        int child2[MAX_CITIES];

        crossover(populationPtrs[j * 2], populationPtrs[j * 2 + 1], child1, cityCount);
        crossover(populationPtrs[j * 2 + 1], populationPtrs[j * 2], child2, cityCount);

        mutate(child1, cityCount);
        mutate(child2, cityCount);

        std::copy(child1, child1 + cityCount, newPopulation[j * 2]);
        std::copy(child2, child2 + cityCount, newPopulation[j * 2 + 1]);
    }

    // Copy the new population back to the main population array
    for (int j = 0; j < MAX_POPULATION; ++j) {
        std::copy(newPopulation[j], newPopulation[j] + cityCount, population[j]);
        populationPtrs[j] = population[j];
    }
}

// After all generations, find the best route
double bestDistance = routeDistance(populationPtrs[0], cityCount);
int bestRoute[MAX_CITIES];
std::copy(populationPtrs[0], populationPtrs[0] + cityCount, bestRoute);

// Gathering the best route from all processes to find the global best
double globalBestDistance;
int globalBestRoute[MAX_CITIES];

// Use custom MPI functions to gather all best distances and routes
 MPI_Allreduce(&bestDistance, &globalBestDistance, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

// Assuming we have a function to collect all routes and find the global best
// This would likely involve gathering all routes and their distances and then selecting the best
my_MPI_Allgather(bestRoute, cityCount, MPI_INT, globalBestRoute, cityCount, MPI_INT, MPI_COMM_WORLD);

if (rank == 0) {
    std::cout << "Global best route distance: " << globalBestDistance << std::endl;
    std::cout << "Global best route: ";
    for (int i = 0; i < cityCount; ++i) {
        std::cout << globalBestRoute[i] << " ";
    }
    std::cout << std::endl;
}

MPI_Finalize();
return 0;
}