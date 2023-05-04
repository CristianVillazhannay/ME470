//Import the header files that we will be using
#include "mpi.h"
#include <math.h>
#include <iostream>
#include <fstream>

//Forward declarations 

//Main Program
int main(int argc, char **argv)
{
	//Serial Code
	
	//Initialize MPI. Begin parallel code.
	MPI_Init(&argc, &argv);

	//Set up the lengths of the domain in all directions.
	double Lx = 1.0;
	double Ly = 1.0;
	double Lz = 1.0;

	//Number of cells in each direction
	int nx = 100;
	int ny = 100;
	int nz = 100;

	//Set up our iteration and tolerance for our code
	int maxIters = 1000; 
	double tolerance = 1.0e-06;

	//Compute the difference in between cells
	double dx = Lx/nx;
	double dy = Ly/ny;
	double dz = Lz/nz;

	//Get my processor ID and number of processors 
	int myProcID, numProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myProcID);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);	

	//Setup our Cartesion topology
	int ndims = 2; 
	int dims[3] = {0, 0, 0};
	MPI_Dims_create(numProcs, ndims, dims);
	int periods[3] = {0, 0, 0};
	int reorder = 2

	//Finalize MPI. End parallel code. 
	MPI_Finalize();

	//Serial Code

	return 0;
}