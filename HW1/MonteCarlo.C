#include "mpi.h"

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>

// Forward declaration of function
double myFunction(double x);
double gen_coord(void);


// Main subroutine
int main(int argc, char **argv)
{
  // Read the total number of points from the command line.
  // The array argv is a size 2 vector because it needs to be able to append the 
  // Program name. This basically assures us that the program will only take one argument. 
  if (argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <numPoints>" << std::endl;
    return 1;
  }

  // Call MPI initialization
  MPI_Init(&argc, &argv);

  // Get my processor ID and number of processor
  int myProcID;
  int numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &myProcID);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // Set number of cells from command line argument. Since argv[] is an array of strings,
  // we need to be able to convert it to an integer in order to be able to use it. 
  int numPoints  = atoi(argv[1]);
  double totalPts = 0; 

  //  Seed the random number generator with each processors unique processor ID.
  srand(myProcID);

  //  Generate the random points

  for (int i = 0; i < numPoints; i++)
  {
    //Random coordinate held
    double x = gen_coord();
    double y = gen_coord();

    //Distance to the centers of both circles
    double cent_dist = sqrt(pow(x-0.5,2) + pow(y-0.5,2));
    double side_dist = sqrt(pow(x-1,2) + pow(y-0.5,2));

    //Needs to be in big circle not small. 
    if (cent_dist <= 0.5 && side_dist >= 0.25)
    {
      totalPts += 1.0;
    } 
  }

  // Sum result across processors
  double totalResult;
  MPI_Reduce(&totalPts, &totalResult, 1, MPI_DOUBLE, MPI_SUM,
             0, MPI_COMM_WORLD);

  // Print the result, but only from root processor
  if (myProcID == 0)
  {
    totalResult = totalResult/ ((double)numProcs * (double)numPoints);
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

double gen_coord(void)
{
  //This will be a normalized version of the RAND_MAX constant, which will generate doubles
  //between 0 and 1.
  return (double)rand() / (double)RAND_MAX;
}

