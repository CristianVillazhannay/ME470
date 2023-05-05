//Import the header files that we will be using
#include "mpi.h"
#include <math.h>
#include <iostream>
#include <fstream>

//Forward declarations 
void getLocalSize(MPI_Comm cartComm, int dims, int n, int &nnLocal, int &offset);
void setToZero(int xSize, int ySize, int zSize, double * u);
void getXYZ(int i, int j, int k, int xOffset, int yOffset, int zOffset,
           double dx, double dy, double dz, double & x, double & y, double & z);
void setBCs(MPI_Comm cartComm, int xSize, int ySize, int zSize, int xOffset, int yOffset,
            int zOffset, double dx, double dy, double dz, double * u);
double jacobiSweep(int xSize, int ySize, int zSize, int xOffset, int yOffset, int zOffset,
                   double dx, double dy, double dz, const double * u, double * unew);
void exchangeBoundaryData(MPI_Comm cartComm, int xSize, int ySize, int zSize, double * u);
void setEqual(int xSize, int ySize, int zSize,const double * u1, double * u2);

// double bcFunc(double x, double y);

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
	

	//Set up the initial state: zeros except for boundary conditions. (Set this for
	//unew too, since the boundary values don't get changed)
  setToZero(gridSizeX, gridSizeY, gridSizeZ, u);
  setToZero(gridSizeX, gridSizeY, gridSizeZ, unew);
  setBCs(cartComm, gridSizeX, gridSizeY, gridSizeZ, xOffset, yOffset, zOffset, dx, dy, dz, u);
  setBCs(cartComm, gridSizeX, gridSizeY, gridSizeZ, xOffset, yOffset, zOffset, dx, dy, dz, unew);

  //Begin the iterations. 
  double dt = 1.0e-4;
  double time = 0.5;
  int maxIters = time / dt;
  int iter = 0; 

  //Iteration Loops
  while (iter < maxIters)
  {
  	//Call the function, and it should update the unew. 
  	jacobiSweep(gridSizeX, gridSizeY, gridSizeZ, xOffset, yOffset, zOffset, 
  							dx, dy, dz, u, unew);

  	//Exchange process data in unew. 
  	exchangeBoundaryData(cartComm, gridSizeX, gridSizeY, gridSizeZ, unew);

  	++iter; 

  	//set u equal to unew.
    setEqual(gridSizeX, gridSizeY, gridSizeZ, unew, u);
  }

  if (myProcID == 0)
  {
  	std::cout << "After " << iter << "iterations, finished" << std::endl;
  }


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
  int cartDims[3];
  int periods[3];
  int cartCoords[3];
  MPI_Cart_get(cartComm, 3, cartDims, periods, cartCoords);

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

void setToZero(int xSize, int ySize, int zSize, double * u)
{
	for (int k = 0; k < zSize; ++k)
	{
	  for (int j = 0; j < ySize; ++j)
	  {
	    for (int i = 0; i < xSize; ++i)
	    {
	      u[(xSize * ySize) * k + xSize * j + i] = 0.0;
	    }
	  }
	}
}

void setBCs(MPI_Comm cartComm, int xSize, int ySize, int zSize, int xOffset, int yOffset,
            int zOffset, double dx, double dy, double dz, double * u)
{
  // Get my process coordinates
  int cartDims[3];
  int periods[3];
  int cartCoords[3];
  MPI_Cart_get(cartComm, 3, cartDims, periods, cartCoords);

  double x, y, z;
  
  // Lower BC --- XZ Plane. 
  if (cartCoords[1] == 0)
  {
  	//Initialize the floor variable (y = 0)
    int j = 0; 

    //Sweep over all the z and the x variables.
    for (int k = 0; k < zSize; ++k)
    {
    	for (int i = 0; i < xSize; ++i)
    	{
    		getXYZ(i, j, k, xOffset, yOffset, zOffset, dy, dy, dz, x, y, z);
    		u[(xSize * ySize) * k + xSize * j + i] = 1.0;
    	}
    }
  }

  // Upper BC -- XZ Plane
  if (cartCoords[1] == cartDims[1] - 1) //Last dimension in the array. 
  {
  	//Initialize the ceiling variable (y = Top)
    int j = ySize - 1; 

    //Sweep over all the z and the x variables.
    for (int k = 0; k < zSize; ++k)
    {
    	for (int i = 0; i < xSize; ++i)
    	{
    		getXYZ(i, j, k, xOffset, yOffset, zOffset, dy, dy, dz, x, y, z);
    		u[(xSize * ySize) * k + xSize * j + i] = 1.0;
    	}
    }
  }

  // Left BC -- YZ Plane
  if (cartCoords[0] == 0)
  {
  	//Intialize the boundary variable (x = 0)
    int i = 0;

    //Sweep over all the z and the x variables.
    for (int k = 0; k < zSize; ++k)
    {
    	for (int j = 0; j < ySize; ++j)
    	{
    		getXYZ(i, j, k, xOffset, yOffset, zOffset, dy, dy, dz, x, y, z);
    		u[(xSize * ySize) * k + xSize * j + i] = 1.0;
    	}
    }
  }

  // Right BC -- YZ Plane. 
  if (cartCoords[0] == cartDims[0] - 1)
  {
  	//Intialize the boundary value (x = Right)
    int i = xSize - 1;

    //Sweep over all the z and the x variables.
    for (int k = 0; k < zSize; ++k)
    {
    	for (int j = 0; j < ySize; ++j)
    	{
    		getXYZ(i, j, k, xOffset, yOffset, zOffset, dy, dy, dz, x, y, z);
    		u[(xSize * ySize) * k + xSize * j + i] = 1.0;
    	}
    }
  }

  //Front BC -- XY Plane
  if (cartCoords[2] == cartDims[2] - 1)
  {
  	//Intialize the boundary value (x = Right)
    int k = 0;

    //Sweep over all the z and the x variables.
    for (int j = 0; j < ySize; ++j)
    {
    	for (int i = 0; i < xSize; ++i)
    	{
    		getXYZ(i, j, k, xOffset, yOffset, zOffset, dy, dy, dz, x, y, z);
    		u[(xSize * ySize) * k + xSize * j + i] = 1.0;
    	}
    }
  }

  //Back BC -- XY Plane
  if (cartCoords[2] == cartDims[2] - 1)
  {
  	//Intialize the boundary value (x = Right)
    int k = zSize - 1;

    //Sweep over all the z and the x variables.
    for (int j = 0; j < ySize; ++j)
    {
    	for (int i = 0; i < xSize; ++i)
    	{
    		getXYZ(i, j, k, xOffset, yOffset, zOffset, dy, dy, dz, x, y, z);
    		u[(xSize * ySize) * k + xSize * j + i] = 1.0;
    	}
    }
  }
}

void getXYZ(int i, int j, int k, int xOffset, int yOffset, int zOffset,
           double dx, double dy, double dz, double & x, double & y, double & z)
{
  x = (i+xOffset) * dx;
  y = (j+yOffset) * dy;
  z = (k+zOffset) * dz;
}

double jacobiSweep(int xSize, int ySize, int zSize, int xOffset, int yOffset, int zOffset,
                   double dx, double dy, double dz, const double * u, double * unew)
{
  double dx2 = dx*dx;
  double dy2 = dy*dy;
  double dz2 = dz*dz;

  double x, y, z;

  for (int k = 1; k < zSize - 1; ++k)
	{
	  for (int j = 1; j < ySize - 1; ++j)
	  {
	    for (int i = 1; i < xSize - 1; ++i)
	    {
	    	double uA = u[(xSize * ySize) * (k    ) + xSize * (j    ) + (i    )];
	    	double uL = u[(xSize * ySize) * (k    ) + xSize * (j    ) + (i - 1)];
	    	double uR = u[(xSize * ySize) * (k    ) + xSize * (j    ) + (i + 1)];
	    	double uB = u[(xSize * ySize) * (k    ) + xSize * (j - 1) + (i    )];
	    	double uT = u[(xSize * ySize) * (k    ) + xSize * (j + 1) + (i    )];
	    	double uF = u[(xSize * ySize) * (k + 1) + xSize * (j    ) + (i    )];
	    	double ub = u[(xSize * ySize) * (k - 1) + xSize * (j    ) + (i    )];

	      unew[(xSize * ySize) * k + xSize * j + i] = 0.1 * ((uR - 2 * uA + uL)/(dx2) + (uT - 2 * uA + uB)/(dy2) + (uF - 2 * uA + ub)/(dz2));
	    }
	  }
	}

  return 0;
}

void exchangeBoundaryData(MPI_Comm cartComm, int xSize, int ySize, int zSize, double * u)
{
  // Define XY, XZ, and ZY Plane datatypes.
  MPI_Datatype XZ, YZ, XY; 

  //This can be contiguous since the MPI can read the data across the 
  //X row then up the Y rows.
  MPI_Type_contiguous(xSize * ySize, MPI_DOUBLE, &XY);

  //Creation of the YZ datatype. This will be a vector type that simply 
  //reads up the column by skipping.
  MPI_Type_vector(xSize * ySize, 1, xSize, MPI_DOUBLE, &YZ);

  //Creation of the XZ datatype. This will be a vector type that skips 
  //up the entirity of the XZ plane to get back to the base row. 
  MPI_Type_vector(xSize * ySize, xSize, xSize * ySize, MPI_DOUBLE, &XZ);

  //Don't forget to commit your changes to the actual variable.
  MPI_Type_commit(&XY);
  MPI_Type_commit(&YZ);
  MPI_Type_commit(&XZ);


  //Get the processes at the front, back, right left, up, and down. 
  int procR, procL, procU, procD, procF, procB;
  MPI_Cart_shift(cartComm, 0, 1, &procL, &procR);
  MPI_Cart_shift(cartComm, 1, 1, &procD, &procR);
  MPI_Cart_shift(cartComm, 2, 1, &procB, &procF);

  MPI_Status status;

  //Communicate up (send my top surface up, recv my bottom surface from down.)
  MPI_Sendrecv(&u[xSize * (ySize -2)], 1, XZ, procU, 0,
  						 &u[0], 1, XZ, procD, 0,
  						 cartComm, &status);

  //Communicate down (send my bottom surface down, recv my top surface from up.)
  MPI_Sendrecv(&u[xSize], 1, XZ, procD, 0,
  						 &u[xSize * (ySize -1)], 1, XZ, procU, 0,
  						 cartComm, &status);

  //Communicate right (send my right surface right, recv my left surface from left.)
  MPI_Sendrecv(&u[xSize - 2], 1, YZ, procR, 0,
  						 &u[0], 1, YZ, procL, 0,
  						 cartComm, &status);

  //Communicate left (send my left surface left, recv my left boundary from right.)
  MPI_Sendrecv(&u[1], 1, YZ, procL, 0,
  						 &u[xSize - 1], 1, YZ, procR, 0,
  						 cartComm, &status);

  //Communicate forward (send my front surface forward, recv my back boundar from behind.)
  MPI_Sendrecv(&u[2*(xSize * ySize)], 1, XY, procF, 0,
  						 &u[(zSize-1)*(xSize*ySize)], 1, XY, procB, 0,
  						 cartComm, &status);

  //Communicate back (send my back surface back, recv my front boundary from forward.)
  MPI_Sendrecv(&u[(zSize-2)*(xSize*ySize)], 1, XZ, procB, 0,
  						 &u[0], 1, XZ, procF, 0,
  						 cartComm, &status);

  //Free all the memory that was done.
  MPI_Type_free(&XZ);
  MPI_Type_free(&YZ);
  MPI_Type_free(&XY); 
}

void setEqual(int xSize, int ySize, int zSize,const double * u1, double * u2)
{
  for (int k = 0; k < zSize; ++k)
	{
	  for (int j = 0; j < ySize; ++j)
	  {
	    for (int i = 0; i < xSize; ++i)
	    {
	      u2[(xSize * ySize) * k + xSize * j + i] = u1[(xSize * ySize) * k + xSize * j + i];
	    }
	  }
	}
}




