#include <petscksp.h>
#include <petscviewer.h>
#include <petscsys.h>

// Forward declarations
double bcFunc(double x, double y);
double forceFunc(double x, double y);

static char help[] = "Use PETSc to solve Poisson's equation in 2D in parallel.\n";

int main(int argc, char **argv)
{
	
  PetscInitialize(&argc, &argv, (char*)0, help);

  // Set some parameters
  double Lx = 1.0;
  double Ly = 1.0;

  int nx = 50; // Number of cells in x direction (nxNodes = nx+1)
  int ny = 50; // Number of cells in y direction (nyNodes = ny+1);

  //int maxIters = 1000;
  double tolerance = 1.0e-6;
  
  // Compute some things
  double dx = Lx/nx;
  double dy = Ly/ny;
  double dx2 = dx * dx;
  double dy2 = dy * dy;

  // The total number of nodes (including boundaries) in each direction:
  int nxNodes = nx + 1;
  int nyNodes = ny + 1;
  int nNodes = nxNodes * nyNodes;
  
  // Create PETSc Vector to hold solution
  Vec u;
  VecCreate(PETSC_COMM_WORLD, &u);
  VecSetSizes(u, PETSC_DECIDE, nNodes);
  VecSetFromOptions(u);

  // Also create PETSc Vector to hold forcing (or boundary conditions)
  // Just duplicate structure of u
  Vec f;
  VecDuplicate(u, &f);

  // Create PETSc A matrix, global size nNodes x nNodes
  Mat A;
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, nNodes, nNodes);
  MatSetFromOptions(A);
  // Preallocate matrix storage (similar to example ex2)
  MatMPIAIJSetPreallocation(A,5,NULL,5,NULL);
  MatSeqAIJSetPreallocation(A,5,NULL);

  // Poisson equation: d2u/dx2 + d2u/dy2 = f(x,y)
  //
  // We are forming a matrix problem A*u = f, where each row of the
  // matrix equation corresponds to the equation at a single node.
  //
  //   For interior nodes, this equation is:
  //     (u_im,j - 2u_i,j + u_ip,j)/dx2 + (u_i,jm - 2u_i,j + u_i,jp)/dy2 = f_i,j
  //       or
  //     1/dx2 * (u_im,j + u_ip,j) + 1/dy2 * (u_i,jm + u_i,jp) - 2*(1/dx2 + 1/dy2)*u_i,j = f_i,j
  //
  //   For boundary nodes, the equation is:
  //     u_i,j = g_i,j where g_i,j is the given boundary value function
  //
  //   Note that this method for applying the BCs destroys matrix symmetry!!!
  //
  // Loop over nodes by looping over the rows of the matrix and assembling.
  // But don't confuse rows of the matrix (1 to nNodes) with rows of the grid (1 to nxNodes)!!!
  //
  // To avoid confusion, let I be the number of the node (matrix row),
  // while i and j are the grid rown and column indices.
  //
  // Only loop over locally owned nodes (obtained from MatGetOwnershipRange)
  //
  PetscInt IStart, IEnd;
  MatGetOwnershipRange(A, &IStart, &IEnd);
  for (int I = IStart; I < IEnd; ++I)
  {
    // Convert node number to i,j indices. The "fast" index is i.
    int i = I % nxNodes;
    int j = I / nxNodes;
    // Get x and y values
    double x = i*dx;
    double y = j*dy;
    // For boundary nodes, u_i,j = g_i,j where g_i,j is the given boundary value function
    // So just put a 1 on the diagonal, set f_i,j = g_i,j
    if ((i == 0) || (i == nxNodes-1) || (j == 0) || (j == nyNodes-1))
    {
      PetscScalar val = 1.0;
      MatSetValues(A, 1, &I, 1, &I, &val, INSERT_VALUES);
      val = bcFunc(x, y);
      VecSetValues(f, 1, &I, &val, INSERT_VALUES);
    }
    // Otherwise, insert proper values for this node's equation:
    //   1/dx2 * (u_im,j + u_ip,j) + 1/dy2 * (u_i,jm + u_i,jp) - 2*(1/dx2 + 1/dy2)*u_i,j = f_i,j
    else
    {
      // Figure out node index for each i-j pair needed
      int Iij  = I;
      int Iimj = I - 1;
      int Iipj = I + 1;
      int Iijm = I - nxNodes;
      int Iijp = I + nxNodes;
      int cols[5] = {Iimj, Iipj, Iijm, Iijp, Iij};
      PetscScalar vals[5] = {1./dx2, 1./dx2, 1./dy2, 1./dy2, -2*(1./dx2 + 1./dy2)};
      MatSetValues(A, 1, &I, 5, cols, vals, INSERT_VALUES);
      PetscScalar val = forceFunc(x, y);
      VecSetValues(f, 1, &I, &val, INSERT_VALUES);
    }
  }

  // Finish assembly
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(f);
  VecAssemblyEnd(f);

  // Create KSP and solve
  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, A, A);
  KSPSetTolerances(ksp, tolerance, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
  KSPSetFromOptions(ksp);
  KSPSolve(ksp, f, u);
  KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);
  
  // Output to text file
  PetscObjectSetName((PetscObject)u, "u");
  PetscViewer viewer;
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "poissonVariables.m", &viewer);
  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  VecView(u, viewer);

  // Clean up
  VecDestroy(&u);
  VecDestroy(&f);
  MatDestroy(&A);
  KSPDestroy(&ksp);
  
  PetscFinalize();
  
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

