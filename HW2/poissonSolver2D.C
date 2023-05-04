#include "mpi.h"
#include <math.h>
#include <iostream>
#include <fstream>

// Forward declarations
void getLocalSize(MPI_Comm cartComm, int dim, int n, int & nLocal, int & offset);
void setToZero(int xSize, int ySize, double * u);
void getXY(int i, int j, int xOffset, int yOffset,
           double dx, double dy, double & x, double & y);
void setBCs(MPI_Comm cartComm, int xSize, int ySize, int xOffset, int yOffset,
            double dx, double dy, double * u);
double bcFunc(double x, double y);
double forceFunc(double x, double y);
double jacobiSweep(int xSize, int ySize, int xOffset, int yOffset,
                   double dx, double dy, const double * u, double * unew);
void setEqual(int xSize, int ySize, const double * u1, double * u2);
void exchangeBoundaryData(MPI_Comm cartComm, int xSize, int ySize, double * u);
void outputToFile(MPI_Comm cartComm, double * u, int nx, int ny, int gridSizeX, int gridSizeY,
                 int xOffset, int yOffset, double dx, double dy);
void writeToFile(const double * u, int nx, int ny, double dx, double dy);


int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  // Set some parameters
  double Lx = 1.0;
  double Ly = 1.0;
  int nx = 50; // Number of cells in x direction (nnodesx = nx+1)
  int ny = 50; // Number of cells in y direction (nnodesy = ny+1);

  int maxIters = 1000;
  double tolerance = 1.0e-6;
  
  // Compute some things
  double dx = Lx/nx;
  double dy = Ly/ny;
  
  // Get my processor ID and number of processors
  int myProcID, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &myProcID);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // Set up Cartesian topology
  int ndims = 2;
  int dims[2] = {0, 0};
  MPI_Dims_create(numProcs, ndims, dims);
  int periods[2] = {0, 0};
  int reorder = 1;
  MPI_Comm cartComm;
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims,
                  periods, reorder, &cartComm);
  int cartRank;
  MPI_Comm_rank(cartComm, &cartRank);

  // Compute size of my local grid.
  int nnxLocal, nnyLocal; // Number of local nodes
  int xOffset, yOffset;
  getLocalSize(cartComm, 0, nx, nnxLocal, xOffset);
  getLocalSize(cartComm, 1, ny, nnyLocal, yOffset);

  // Allocate memory for grid -- (nnxLocal+2) by (nnyLocal+2) to account
  // for boundary or ghost nodes  
  int gridSizeX = nnxLocal + 2;
  int gridSizeY = nnyLocal + 2;
  double * u    = new double[gridSizeX * gridSizeY];
  double * unew = new double[gridSizeX * gridSizeY];
  
  // Set up initial state: zeros except for boundary condition. (Set
  // this for unew, too, since the boundary values don't get changed
  // durint the Jacobi sweep step.)
  setToZero(gridSizeX, gridSizeY, u);
  setToZero(gridSizeX, gridSizeY, unew);
  setBCs(cartComm, gridSizeX, gridSizeY, xOffset, yOffset, dx, dy, u);
  setBCs(cartComm, gridSizeX, gridSizeY, xOffset, yOffset, dx, dy, unew);
  
  // Begin Jacobi loop
  int iter = 0;
  double diff = 2*tolerance;
  while (iter < maxIters && diff > tolerance)
  {

    // Jacobi sweep
    double localDiff = jacobiSweep(gridSizeX, gridSizeY, xOffset, yOffset,
                                   dx, dy, u, unew);

    // Exchange process boundary data in unew
    exchangeBoundaryData(cartComm, gridSizeX, gridSizeY, unew);
    
    // Check for convergence
    // "localDiff" returned the sum of the squared differences over the nodes.
    // Reduce, take the square root, and normalize by number of interior nodes.
    MPI_Allreduce(&localDiff, &diff, 1, MPI_DOUBLE, MPI_SUM, cartComm);
    diff = sqrt(diff)/((nx-1) * (ny-1));
    ++iter;

    // Set u equal to unew
    setEqual(gridSizeX, gridSizeY, unew, u);

  }

  if (myProcID == 0)
  {
    std::cout << "After " << iter << " iterations, finished with du = " << diff << std::endl;
  }

  // Output
  outputToFile(cartComm, u, nx, ny, gridSizeX, gridSizeY, xOffset, yOffset, dx, dy);
  
  // Clean up
  delete[] u;
  delete[] unew;
  
  
  MPI_Finalize();
  return 0;
}

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

void setToZero(int xSize, int ySize, double * u)
{
  for (int j = 0; j < ySize; ++j)
  {
    for (int i = 0; i < xSize; ++i)
    {
      u[xSize * j + i] = 0.0;
    }
  }
}

void setBCs(MPI_Comm cartComm, int xSize, int ySize, int xOffset, int yOffset,
            double dx, double dy, double * u)
{
  // Get my process coordinates
  int cartDims[2];
  int periods[2];
  int cartCoords[2];
  MPI_Cart_get(cartComm, 2, cartDims, periods, cartCoords);

  double x, y;
  
  // Lower BC?
  if (cartCoords[1] == 0)
  {
    int j = 0; 
    for (int i = 0; i < xSize; ++i)
    {
      getXY(i, j, xOffset, yOffset, dx, dy, x, y);
      u[xSize * j + i] = bcFunc(x,y);
    }
  }

  // Upper BC?
  if (cartCoords[1] == cartDims[1] - 1)
  {
    int j = ySize - 1; 
    for (int i = 0; i < xSize; ++i)
    {
      getXY(i, j, xOffset, yOffset, dx, dy, x, y);
      u[xSize * j + i] = bcFunc(x,y);
    }
  }

  // Left BC?
  if (cartCoords[0] == 0)
  {
    int i = 0;
    for (int j = 0; j < ySize; ++j)
    {
      getXY(i, j, xOffset, yOffset, dx, dy, x, y);
      u[xSize * j + i] = bcFunc(x,y);
    }
  }

  // Right BC?
  if (cartCoords[0] == cartDims[0] - 1)
  {
    int i = xSize - 1;
    for (int j = 0; j < ySize; ++j)
    {
      getXY(i, j, xOffset, yOffset, dx, dy, x, y);
      u[xSize * j + i] = bcFunc(x,y);
    }
  }

}

void getXY(int i, int j, int xOffset, int yOffset,
           double dx, double dy, double & x, double & y)
{
  x = (i+xOffset) * dx;
  y = (j+yOffset) * dy;
}

double bcFunc(double x, double y)
{
  const double pi = 4*atan(1.0);
  double u = cos(pi*x) * sin(3*pi*y);
  return u;
}

double forceFunc(double x, double y)
{
  const double pi = 4*atan(1.0);
  double f = -10 * pi * pi * cos(pi*x) * sin(3*pi*y);
  return f;
}

double jacobiSweep(int xSize, int ySize, int xOffset, int yOffset,
                   double dx, double dy, const double * u, double * unew)
{
  double dx2 = dx*dx;
  double dy2 = dy*dy;
  double x, y;
  double diff2Total = 0.0;
  for (int j = 1; j < ySize - 1; ++ j)
  {
    for (int i = 1; i < xSize - 1; ++i)
    {
      double uE = u[xSize * (j  ) + i+1 ];
      double uW = u[xSize * (j  ) + i-1 ];
      double uN = u[xSize * (j+1) + i   ];
      double uS = u[xSize * (j-1) + i   ];
      getXY(i, j, xOffset, yOffset, dx, dy, x, y);
      double f = forceFunc(x,y);

      unew[xSize * j + i] = (dy2*(uE+uW) + dx2*(uN+uS) - dx2*dy2*f)/(2 * (dx2 + dy2));

      double du = unew[xSize * j + i] - u[xSize * j + i];
      diff2Total += du*du;
    }
  }
  return diff2Total;
}

void setEqual(int xSize, int ySize, const double * u1, double * u2)
{
  for (int j = 0; j < ySize; ++ j)
  {
    for (int i = 0; i < xSize; ++i)
    {
      u2[xSize * j + i] = u1[xSize * j + i];
    }
  }
}

void exchangeBoundaryData(MPI_Comm cartComm, int xSize, int ySize, double * u)
{

  // Define row and column data types
  MPI_Datatype rowType, colType;
  MPI_Type_contiguous(xSize, MPI_DOUBLE, &rowType);
  MPI_Type_vector(ySize, 1, xSize, MPI_DOUBLE, &colType);
  MPI_Type_commit(&rowType);
  MPI_Type_commit(&colType);

  // Get the processes at right, left, up and down
  int procR, procL, procU, procD;
  MPI_Cart_shift(cartComm, 0, 1, &procL, &procR);
  MPI_Cart_shift(cartComm, 1, 1, &procD, &procU);

  MPI_Status status;
  
  // Communicate up (send my top row up, recv my bottom boundary from down)
  MPI_Sendrecv(&u[xSize * (ySize-2)], 1, rowType, procU, 0,
               &u[0],                 1, rowType, procD, 0,
               cartComm, &status);

  // Communicate down (send my bottom row down, recv my top boundary from up)
  MPI_Sendrecv(&u[xSize],             1, rowType, procD, 0,
               &u[xSize * (ySize-1)], 1, rowType, procU, 0,
               cartComm, &status);

  // Communicate right (send my right col right, recv my left boundary from left)
  MPI_Sendrecv(&u[xSize-2], 1, colType, procR, 0,
               &u[0],       1, colType, procL, 0,
               cartComm, &status);
  
  // Communicate left (send my left col left, recv my right boundary from right)
  MPI_Sendrecv(&u[1],       1, colType, procL, 0,
               &u[xSize-1], 1, colType, procR, 0,
               cartComm, &status);

  MPI_Type_free(&rowType);
  MPI_Type_free(&colType);
  
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
