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
double jacobiSweep(double dt, int xSize, int ySize, int zSize, int xOffset, int yOffset, int zOffset,
                   double dx, double dy, double dz, const double * u, double * unew);
void exchangeBoundaryData(MPI_Comm cartComm, int xSize, int ySize, int zSize, double * u);
void setEqual(int xSize, int ySize, int zSize,const double * u1, double * u2);
void outputToFile(MPI_Comm cartComm, double * u, int nx, int ny, int gridSizeX, int gridSizeY,
                 int xOffset, int yOffset, double dx, double dy);
void writeToFile(const double * u, int nx, int ny, double dx, double dy);
void dataWrite(double u);
// void outputToFile(MPI_Comm cartComm, double * u, int nx, int ny, int nz, int gridSizeX, int gridSizeY,
//                   int gridSizeZ, int xOffset, int yOffset, int zOffset, double dx, double dy, double dz);
// void writeToFile(const double * u, int nx, int ny, int nz, double dx, double dy, double dz);


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

	//Compute the difference in between cells
	double dx = Lx/nx;
	double dy = Ly/ny;
	double dz = Lz/nz;

	//Get my processor ID and number of processors 
	int myProcID, numProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myProcID);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);	

	//Setup our Cartesian topology
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
  	jacobiSweep(dt, gridSizeX, gridSizeY, gridSizeZ, xOffset, yOffset, zOffset, 
  							dx, dy, dz, u, unew);

  	//Exchange process data in unew. 
  	exchangeBoundaryData(cartComm, gridSizeX, gridSizeY, gridSizeZ, unew);

  	++iter; 

  	//set u equal to unew.
    setEqual(gridSizeX, gridSizeY, gridSizeZ, unew, u);
  }

  //outputToFile(cartComm, u_slice, nx, ny, gridSizeX, gridSizeY, xOffset, yOffset, dx, dy);
  //outputToFile(cartComm, u, nx, ny, nz, gridSizeX, gridSizeY, gridSizeZ, xOffset, yOffset, zOffset, dx, dy, dz);

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
  	//Intialize the boundary value (z = 0)
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
  	//Intialize the boundary value (z = back)
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

double jacobiSweep(double dt, int xSize, int ySize, int zSize, int xOffset, int yOffset, int zOffset,
                   double dx, double dy, double dz, const double * u, double * unew)
{
  double dx2 = dx*dx;
  double dy2 = dy*dy;
  double dz2 = dz*dz;

  double x, y, z;
  
  //we avoid the edges of the grid.
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

	      unew[(xSize * ySize) * k + xSize * j + i] = uA + (double)0.1 * dt *((uR - (double)2 * uA + uL)/(dx2) + (uT - (double)2 * uA + uB)/(dy2) + (uF - (double)2 * uA + ub)/(dz2));

	      getXYZ(i, j, k, xOffset, yOffset, zOffset, dx, dy, dz, x, y, z);

        if (x == 0.5 && y == 0.5 && z == 0.5) {
          //std::cout << "uA: " << uA << std::endl;
          dataWrite(uA);
        }

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

  //Note that for these datatypes, the count * blocklength should be equal to the
  //number of elements in the desired plane.

  //Creation of the YZ datatype. This will be a vector type that simply 
  //reads up the column by skipping up through the xSize columns.
  MPI_Type_vector(zSize * ySize,     1, xSize, MPI_DOUBLE, &YZ);

  //Creation of the XZ datatype. This will be a vector type that skips 
  //up the entirity of the XZ plane to get back to the base row. 
  MPI_Type_vector(zSize, xSize, xSize * ySize, MPI_DOUBLE, &XZ);

  //Don't forget to commit your changes to the actual variable.
  MPI_Type_commit(&XY);
  MPI_Type_commit(&YZ);
  MPI_Type_commit(&XZ);


  //Get the processes at the front, back, right left, up, and down. This will only 
  //tell us the rank of the processes in the topology.
  int procR, procL, procU, procD, procF, procB;
  MPI_Cart_shift(cartComm, 0, 1, &procL, &procR);
  MPI_Cart_shift(cartComm, 1, 1, &procD, &procU);
  MPI_Cart_shift(cartComm, 2, 1, &procB, &procF);

  MPI_Status status;

  //Communicate up (send my top surface up, recv my bottom surface from down.)
  MPI_Sendrecv(&u[xSize * (ySize - 2)], 1, XZ, procU, 0,
  						 &u[0], 								  1, XZ, procD, 0,
  						 cartComm, &status);

  //Communicate down (send my bottom surface down, recv my top surface from up.)
  MPI_Sendrecv(&u[xSize], 						 1, XZ, procD, 0,
  						 &u[xSize * (ySize -1)], 1, XZ, procU, 0,
  						 cartComm, &status);

  //Communicate right (send my right surface right, recv my left surface from left.)
  MPI_Sendrecv(&u[xSize - 2], 1, YZ, procR, 0,
  						 &u[0], 				1, YZ, procL, 0,
  						 cartComm, &status);

  //Communicate left (send my left surface left, recv my left boundary from right.)
  MPI_Sendrecv(&u[1], 				1, YZ, procL, 0,
  						 &u[xSize - 1], 1, YZ, procR, 0,
  						 cartComm, &status);

  //Communicate forward (send my front surface forward, recv my back boundary from behind.)
  MPI_Sendrecv(&u[xSize * ySize], 					1, XY, procF, 0,
  						 &u[(zSize-1)*(xSize*ySize)], 1, XY, procB, 0,
  						 cartComm, &status);

  //Communicate back (send my back surface back, recv my front boundary from forward.)
  MPI_Sendrecv(&u[(zSize-2)*(xSize*ySize)], 1, XY, procB, 0,
  						 &u[0], 											1, XY, procF, 0,
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

void outputToFile(MPI_Comm cartComm, double * u, int nx, int ny, int gridSizeX, int gridSizeY,
                  int xOffset, int yOffset, double dx, double dy)
{
  int myProcID, numProcs;
  MPI_Comm_rank(cartComm, &myProcID);
  MPI_Comm_size(cartComm, &numProcs);

  // On process 0, allocate data for entire u array, as well as arrays
  // to hold grid sizes and offsets gathered from individual procs
  double * uAll;
  int *gridSizeXArray, *gridSizeYArray, *xOffsetArray, *yOffsetArray;
  int fullSizeX = nx + 1;
  int fullSizeY = ny + 1;
  if (myProcID == 0)
  {
    uAll = new double[fullSizeX*fullSizeY];
    gridSizeXArray = new int[numProcs];
    gridSizeYArray = new int[numProcs];
    xOffsetArray = new int[numProcs];
    yOffsetArray = new int[numProcs];
  }

  // Gather grid sizes and offsets
  MPI_Gather(&gridSizeX, 1, MPI_INTEGER,
             gridSizeXArray, 1, MPI_INTEGER, 0, cartComm);
  MPI_Gather(&gridSizeY, 1, MPI_INTEGER,
             gridSizeYArray, 1, MPI_INTEGER, 0, cartComm);

  MPI_Gather(&xOffset, 1, MPI_INTEGER,
             xOffsetArray, 1, MPI_INTEGER, 0, cartComm);
  MPI_Gather(&yOffset, 1, MPI_INTEGER,
             yOffsetArray, 1, MPI_INTEGER, 0, cartComm);

  // On each processor, send grid data to process 0. Use non-blocking
  // communication to avoid deadlock.
  MPI_Request request;
  MPI_Isend(u, gridSizeX*gridSizeY, MPI_DOUBLE, 0, 0, cartComm, &request);

  // On process 0, loop over processes and receive the sub-block using
  // a derived data type
  if (myProcID == 0)
  {
    for (int proc = 0; proc < numProcs; ++proc)
    {
      MPI_Datatype subblockType;
      int count = gridSizeYArray[proc];
      int length = gridSizeXArray[proc];
      int stride = fullSizeX;
      MPI_Type_vector(count, length, stride, MPI_DOUBLE, &subblockType);
      MPI_Type_commit(&subblockType);
      double * recvPointer = &uAll[ yOffsetArray[proc] * fullSizeX + xOffsetArray[proc] ];
      MPI_Status status;
      MPI_Recv(recvPointer, 1, subblockType, proc,
               0, cartComm, &status);
      MPI_Type_free(&subblockType);
    }
  }

  MPI_Wait(&request, MPI_STATUS_IGNORE);
  
  // Output to file from proc 0
  if (myProcID == 0)
  {
    writeToFile(uAll, fullSizeX, fullSizeY, dx, dy);
  }
  
  // Delete arrays on proc 0
  if (myProcID == 0)
  {
    delete[] uAll;
    delete[] gridSizeXArray;
    delete[] gridSizeYArray;
    delete[] xOffsetArray;
    delete[] yOffsetArray;
  }
  
}

void writeToFile(const double * u, int nx, int ny, double dx, double dy)
{
  std::ofstream myFile;
  myFile.open("poisson2D.dat");

  const int width = 8;
  const int prec = 5;

  for (int j = 0; j < ny; ++j)
  {
    double y = j*dy;
    for (int i = 0; i < nx; ++i)
    {
      double x = i*dx;
      double val = u[j*nx + i];
      char buf[256];
      char pattern[] = "%8.5f\t%8.5f\t%8.5f";
      sprintf(buf, pattern, x, y, val);
      myFile << buf << std::endl;
    }
  }

  myFile.close();
  
}

void dataWrite(double u)
{
  std::ofstream myFile;
  //append not overwrite. This means that middle.dat needs to be deleted everytime this is run.
  myFile.open("middle.dat", std::ios_base::app);

  //Change the precision to be able to read at least 6 sig. figs.
  char buf[256];
  char pattern[] = "%1.8f\t"; 
  sprintf(buf, pattern, u);
  myFile << buf << std::endl;

  myFile.close();
  
}

// void outputToFile(MPI_Comm cartComm, double * u, int nx, int ny, int nz, int gridSizeX, int gridSizeY,
//                   int gridSizeZ, int xOffset, int yOffset, int zOffset, double dx, double dy, double dz)
// {
//   int myProcID, numProcs;
//   MPI_Comm_rank(cartComm, &myProcID);
//   MPI_Comm_size(cartComm, &numProcs);

//   // On process 0, allocate data for entire u array, as well as arrays
//   // to hold grid sizes and offsets gathered from individual procs
//   double * uAll;
//   int *gridSizeXArray, *gridSizeYArray, *gridSizeZArray, *xOffsetArray, *yOffsetArray, *zOffsetArray;
//   int fullSizeX = nx + 1;
//   int fullSizeY = ny + 1;
//   int fullSizeZ = nz + 1;
//   if (myProcID == 0)
//   {
//     uAll = new double[fullSizeX*fullSizeY*fullSizeZ];
//     gridSizeXArray = new int[numProcs];
//     gridSizeYArray = new int[numProcs];
//     gridSizeZArray = new int[numProcs];

//     xOffsetArray = new int[numProcs];
//     yOffsetArray = new int[numProcs];
//     zOffsetArray = new int[numProcs];
//   }

//   // Gather grid sizes and offsets
//   MPI_Gather(&gridSizeX, 1, MPI_INTEGER,
//              gridSizeXArray, 1, MPI_INTEGER, 0, cartComm);
//   MPI_Gather(&gridSizeY, 1, MPI_INTEGER,
//              gridSizeYArray, 1, MPI_INTEGER, 0, cartComm);
//   MPI_Gather(&gridSizeZ, 1, MPI_INTEGER,
//   					 gridSizeZArray, 1, MPI_INTEGER, 0, cartComm);


//   MPI_Gather(&xOffset, 1, MPI_INTEGER,
//              xOffsetArray, 1, MPI_INTEGER, 0, cartComm);
//   MPI_Gather(&yOffset, 1, MPI_INTEGER,
//              yOffsetArray, 1, MPI_INTEGER, 0, cartComm);
//   MPI_Gather(&zOffset, 1, MPI_INTEGER,
//   					 zOffsetArray, 1, MPI_INTEGER, 0, cartComm);

//   // On each processor, send grid data to process 0. Use non-blocking
//   // communication to avoid deadlock.
//   MPI_Request request;
//   MPI_Isend(u, gridSizeX*gridSizeY*gridSizeZ, MPI_DOUBLE, 0, 0, cartComm, &request);

//   // On process 0, loop over processes and receive the sub-block using
//   // a derived data type
//   if (myProcID == 0)
//   {
//     for (int proc = 0; proc < numProcs; ++proc)
//     {
//     	//This sub-block type is made specifically to handle the 2D portion of the code
//     	//where we move up the face of the actual data domain. 
//       MPI_Datatype subblockType;
//       int count = gridSizeYArray[proc];
//       int length = gridSizeXArray[proc];
//       int stride = fullSizeX;
//       MPI_Type_vector(count, length, stride, MPI_DOUBLE, &subblockType);
//       MPI_Type_commit(&subblockType);

//       //This sub-block type will be made specifically to capture a single row of the data;
//       MPI_Datatype rowType;
//       int rowCount = 1;
//       int rowLength = gridSizeXArray[proc];
//       int rowStride = fullSizeX;
//       MPI_Type_vector(rowCount, rowLength, rowStride, MPI_DOUBLE, &rowType);
//       MPI_Type_commit(&rowType);

//       //This will handle the 3D portion of the domain. Since
//       //the data no longer moves up in a linear fashion due the 3D topology, we 
//       //need to account for the backwards movement when a single plane is split.
//       MPI_Datatype dimType;
//       int Dcount  = gridSizeZArray[proc]; //We need Z blocks since each block will planar.
//       int Dlength = gridSizeYArray[proc]; //Creates the face like in subblockType. 
//       int Dstride = fullSizeY;  					//Stride is the entire plane of rows. 
//       MPI_Type_vector(Dcount, Dlength, Dstride, rowType, &dimType);
//       MPI_Type_commit(&dimType);


//       //Pointer to where wwe will begin to assign the "subblock" data type into the actual array. 
//       double * recvPointer = &uAll[zOffsetArray[proc] * fullSizeY * fullSizeX + yOffsetArray[proc] * fullSizeX + xOffsetArray[proc]];
//       MPI_Status status;

//       MPI_Recv(recvPointer, gridSizeZArray[proc] * gridSizeYArray[proc], rowType, proc,
//                0, cartComm, &status);

//       //Free all committed variables.
//       MPI_Type_free(&subblockType);
//       MPI_Type_free(&rowType);
//       MPI_Type_free(&dimType);
//     }
//   }

//   MPI_Wait(&request, MPI_STATUS_IGNORE);
  
//   // Output to file from proc 0
//   if (myProcID == 0)
//   {
//     writeToFile(uAll, fullSizeX, fullSizeY, fullSizeZ, dx, dy, dz);
//   }
  
//   // Delete arrays on proc 0
//   if (myProcID == 0)
//   {
//     delete[] uAll;
//     delete[] gridSizeXArray;
//     delete[] gridSizeYArray;
//     delete[] xOffsetArray;
//     delete[] yOffsetArray;
//   }
  
// }

// void writeToFile(const double * u, int nx, int ny, int nz, double dx, double dy, double dz)
// {
//   std::ofstream myFile;
//   myFile.open("poisson3D.dat");

//   const int width = 8;
//   const int prec = 5;

//   for (int k = 0; k < nz; ++k)
//   {
//   	double z = k*dz;
// 	  for (int j = 0; j < ny; ++j)
// 	  {
// 	    double y = j*dy;
// 	    for (int i = 0; i < nx; ++i)
// 	    {
// 	      double x = i*dx;
// 	      double val = u[k*(nx*ny)+j*nx + i];
// 	      char buf[256];
// 	      char pattern[] = "%8.5f\t%8.5f\t%8.5f\t%8.5f";
// 	      sprintf(buf, pattern, x, y, z, val);
// 	      myFile << buf << std::endl;
// 	    }
// 	  }
// 	}

//   myFile.close();
  
// }
