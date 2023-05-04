//Import the header files that we will be using
#include "mpi.h"
#include <math.h>
#include <iostream>
#include <fstream>

//Forward declarations 
void getLocalSize(MPI_Comm cartComm, int dims, int n, int &nnLocal, int &offset);

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
	int ndims = 3; 
	int dims[3] = {0, 0, 0};
	MPI_Dims_create(numProcs, ndims, dims);


	int periods[3] = {0, 0, 0}; //The problem is not periodic in any direction.
	int reorder = 1;			//May take advantage of physical architecture.
	MPI_Comm cartComm;
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, 
					periods, reorder, &cartComm);

	//Get the processors' new rank in the communicator. 
	int cartRank;
	MPI_Comm_rank(cartComm, &cartRank);

	//Compute the size of the local grid. Domain decomposition.
	int nnxLocal, nnyLocal, nnzLocal;
	int xOffset, yOffset, zOffset; 
	getLocalSize(cartComm, 0, nx, nnxLocal, xOffset);
	getLocalSize(cartComm, 1, ny, nnyLocal, yOffset);
	getLocalSize(cartComm, 2, nz, nnzLocal, zOffset);	

	//Allocate memory for grid -- (nnxLocal + 2) by (nnyLocal+2) to account for 
	//boundary or ghost nodes.
	int gridSizeX = nnxLocal + 2; 
	int gridSizeY = nnyLocal + 2; 	
	int gridSizeZ = nnzLocal + 2; 
	double * u 		 = new double[gridSizeX * gridSizeY * gridSizeZ];
	double * unew  = new double[gridSizeX * gridSizeY * gridSizeZ];	
	

	//Set up the initial state: zeros except for boundary conditions. ()


	//Deallocate the dynamic memory for the u and unew arrays.
	delete[] u;
	delete[] unew;

	//Finalize MPI. End parallel code. 
	MPI_Finalize();

	//Serial Code

	return 0;
}


//Function Declarations
void getLocalSize(MPI_Comm cartComm, int dim, int n, int & nnLocal, int & offset)
{
  // Get my process coordinates
  int cartDims[2];
  int periods[2];
  int cartCoords[2];
  MPI_Cart_get(cartComm, 2, cartDims, periods, cartCoords);

  // Compute number of local interior nodes, and offset
  // global_I = local_i + offset
  nnLocal = (n-1)/cartDims[dim];
  offset = cartCoords[dim] * nnLocal;
  // Account for non-uniform distribution on last proc: offset + nnlocal + 1 = n
  if (cartCoords[dim] == cartDims[dim] - 1)
  {
    nnLocal = n - 1 - offset;
  }
}

