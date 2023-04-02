#include <stdlib.h>
#include <iostream>
#include <iomanip>

// Forward declaration of function
double myFunction(double x);

// Main subroutine
int main(int argc, char **argv)
{
  // Check to make sure number of cells is set by command line
  if (argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <nCells>" << std::endl;
    return 1;
  }

  // Set number of cells from command line argument
  int nCells  = atoi(argv[1]);

  // Define some variables
  double xMin = 0.0;                   
  double xMax = 1.0;
  double dx   = (xMax - xMin)/nCells;

  // Loop over cells and compute integral using trapezoidal rule
  double result = 0.0;
  for (int i = 0; i < nCells; ++i)
  {
    double xL = xMin + i*dx;
    double xR = xL + dx;
    result += 0.5 * (myFunction(xL)+myFunction(xR)) * dx;
  }

  // Print the result
  std::cout << "result = ";
  std::cout << std::fixed << std::setprecision(10) << result << std::endl;
  
  return 0;
}

// Define function
double myFunction(double x)
{
  double f = 4.0/(1.0 + x*x);
  return f;
}

