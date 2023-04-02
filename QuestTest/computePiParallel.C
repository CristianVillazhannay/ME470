#include "mpi.h"

#include <iostream>
#include <iomanip>
#include <stdlib.h>

// Forward declaration of function
double myFunction(double x);

// Main subroutine
int main(int argc, char **argv)
{
  // Check to make sure number of cells is set by command line
  if (argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <numCells>" << std::endl;
    return 1;
  }

  // Call MPI initialization
  MPI_Init(&argc, &argv);

  // Get my processor ID and number of processor
  int myProcID;
  int numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &myProcID);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  
  // Set number of cells from command line argument
  int numCells  = atoi(argv[1]);

  // Define some variables
  double xMin = 0.0;                   
  double xMax = 1.0;
  double dx   = (xMax - xMin)/numCells;

  // Find which piece of the integral this processor is responsible for
  int numCellsPerProc = numCells/numProcs;
  int myStart = myProcID * numCellsPerProc;
  int myEnd   = myStart + numCellsPerProc;

  // Account for unequal load by making sure last processor has correct end
  if (myProcID == numProcs - 1) myEnd = numCells;
  
  // Loop over cells and compute integral using trapezoidal rule
  double myResult = 0.0;
  for (int i = myStart; i < myEnd; ++i)
  {
    double xL = xMin + i*dx;
    double xR = xL + dx;
    myResult += 0.5 * (myFunction(xL)+myFunction(xR)) * dx;
  }

  // Sum result across processors
  double totalResult;
  MPI_Reduce(&myResult, &totalResult, 1, MPI_DOUBLE, MPI_SUM,
             0, MPI_COMM_WORLD);
  
  // Print the result, but only from root processor
  if (myProcID == 0)
  {
    std::cout << "result = ";
    std::cout << std::fixed << std::setprecision(10)
              << totalResult << std::endl;
  }

  // Call MPI_Finalize
  MPI_Finalize();
  
  return 0;
}

// Define function
double myFunction(double x)
{
  double f = 4.0/(1.0 + x*x);
  return f;
}

