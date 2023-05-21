#include <petscksp.h>
#include <fstream>

static char help[] = "Use PETSc to solve the 3D time-dependent heat equation.\n";

// Forward declarations
void writeValues(int n, const double *t, const PetscScalar *u, const char file[]);

int main(int argc, char **argv)
{
  PetscInitialize(&argc, &argv, (char*)0, help);

  //Initialize the name flag for the output .dat 
  int max = 20;
  char name[20];
  PetscBool set; 
  PetscOptionsGetString(NULL, "-name", name, max, &set);

  // Set some parameters
  double D = 0.1;

  double Lx = 1.0;
  double Ly = 1.0;
  double Lz = 1.0;

  int nx = 100; // Number of cells in x direction (nxNodes = nx+1)
  PetscOptionsGetInt(NULL, "-nx", &nx, NULL);
  int ny = nx; // Number of cells in y direction (nyNodes = ny+1);
  int nz = nx; // Number of cells in z direction (nzNodes = nz+1);

  double dt = 0.01;
  PetscOptionsGetReal(NULL, "-dt", &dt, NULL);
  double totalTime = 0.50;
  
  double tol = 1.0e-6;

  // Compute some things
  double dx = Lx/nx;
  double dy = Ly/ny;
  double dz = Lz/nz;
  double dx2 = dx * dx;
  double dy2 = dy * dy;
  double dz2 = dz * dz;
  int numTimeSteps = totalTime/dt;

  // Total number of nodes (including boundaries) in each direction
  int nxNodes = nx + 1;
  int nyNodes = ny + 1;
  int nzNodes = nz + 1;
  int nNodes = nxNodes * nyNodes * nzNodes;

  // Start timer
  double startTime = MPI_Wtime();
  
  // Create PETSc Vectors to hold solution at timesteps n and n+1, and the increment
  Vec un, unp1, du;
  VecCreate(PETSC_COMM_WORLD, &un);
  VecSetSizes(un, PETSC_DECIDE, nNodes);
  VecSetFromOptions(un);
  VecDuplicate(un, &unp1);
  VecDuplicate(un, &du);

  // Create PETSc Vector to hold RHS of equation; same structure as u
  Vec b;
  VecDuplicate(un, &b);

  // Set un initially to zero (BC will be applied below)
  VecSet(un, 0.0);
  
  // Create PETSc Matrix J, global size nNodes x nNodes
  Mat J;
  MatCreate(PETSC_COMM_WORLD, &J);
  MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, nNodes, nNodes);
  MatSetFromOptions(J);
  // Preallocate matrix storage (similar to example ex2)
  MatMPIAIJSetPreallocation(J,7,NULL,7,NULL);
  MatSeqAIJSetPreallocation(J,7,NULL);

  // Heat equation: 
  //  du/dt = D*(d2u/dx2 + d2u/dy2 + d2u/dz2)
  //
  // Backward Euler time integration:
  //  -- discretize du/dt = (u^(n+1) - u^n)/dt
  //  -- discretize spatial derivatives using n+1 timestep quantities
  // 
  // F_ijk(unp1) = 0 = ?
  //
  //
  // Iterative scheme:
  //   unp1kp1 = unp1k + du
  //
  // Write equation for du in the form:
  //
  // J(unp1k) * du  = -F_ijk(unp1k)
  //
  // Note that if the LHS can be written as J*du, then F_ijk(unp1) = J*unp1 - un
  // and the RHS vector is -F_ijk(unp1) = un - J*unp1
  // We can use this to quickly update the RHS at each step.
  //
  // On boundary nodes, the equation is simply
  //   u_ijk = g_ijk where g_ijk is the given boundary value function (1)
  //
  // Loop over nodes  by looping over the rows of the matrix and assembling.
  // But don't confuse rows of the matrix (1 to nNodes) with rows of the grid (1 to nxNodes)!!!
  //
  // To avoid confusion, let I be the number of the node (matrix row),
  // while i,j,k are the grid indices.
  //
  // Only loop over locally owned nodes (obtained from MatGetOwnershipRange)
  //
  PetscInt IStart, IEnd;
  MatGetOwnershipRange(J, &IStart, &IEnd);
  for (int I = IStart; I < IEnd; ++I)
  {
    // Convert node number to i,j,k indices. The "fast" index is i.
    int k = I / (nxNodes * nyNodes);
    int j = (I - k*nxNodes * nyNodes) / nxNodes;
    int i = I - k * nxNodes * nyNodes - j * nxNodes;
    // For boundary nodes, u_ijk = g_ijk = 1.0
    // So, du_ijk = 0 at each step (no change in u_ijk)
    // Just put a 1 on the diagonal, and also set BC on un
    if ((i == 0) || (i == nxNodes-1) || (j == 0) || (j == nyNodes-1) || (k == 0) || (k == nzNodes-1))
    {
      PetscScalar val = 1.0;
      MatSetValues(J, 1, &I, 1, &I, &val, INSERT_VALUES);
      val = 1.0;
      VecSetValues(un, 1, &I, &val, INSERT_VALUES);
    }    
    // Otherwise, insert coefficients for J matrix
    else
    {
      // Figure out node index for each ijk triplet needed
      int Iijk  = I;
      int Iimjk = I - 1;
      int Iipjk = I + 1;
      int Iijmk = I - nxNodes;
      int Iijpk = I + nxNodes;
      int Iijkm = I - nxNodes * nyNodes;
      int Iijkp = I + nxNodes * nyNodes;

      // TODO: Uncomment and fix the next 3 lines
      int cols[7] = {Iimjk, Iipjk, Iijmk, Iijpk,Iijkm,Iijkp,Iijk};
      PetscScalar vals[7] = {1./dx2, 1./dx2, 1./dy2, 1./dy2, 1./dz2, 1./dz2, -2*(1./dx2 + 1./dy2 + 1./dz2)};
      MatSetValues(J, 1, &I, 7, cols, vals, INSERT_VALUES);
    }
  }
  
  // Finish assembly
  MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(un);
  VecAssemblyEnd(un);

  // Create KSP
  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, J, J);
  KSPSetTolerances(ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
  KSPSetFromOptions(ksp);
  
  // Figure out which processor holds the node at x=y=z=0.5, and allocate storage
  int imid = nx/2;
  int jmid = ny/2;
  int kmid = nz/2;
  int Imid = kmid*nxNodes*nyNodes + jmid*nxNodes + imid;
  bool haveMiddle = (IStart <= Imid) && (IEnd > Imid);
  PetscScalar * uMidValues = NULL;
  double      * timeValues = NULL;
  if (haveMiddle)
  {
    uMidValues = new PetscScalar[numTimeSteps+1];
    timeValues = new double[numTimeSteps+1];
    uMidValues[0] = 0.0;
    timeValues[0] = 0.0;
  }
  
  // Begin time integration loop
  for (int iTime = 0; iTime < numTimeSteps; ++iTime)
  {
    // Set unp1 = un initially
    // TODO: Set unp1 = un
    VecCopy(un, unp1);
    
    // Form RHS Vector
    // This is just b = un - J*unp1
    // (note that this also gives the correct value on the boundary nodes, b = 0
    // TODO: Use PETSc matrix-vector operations to set b = un - J*unp1
    // (See p. 44 and p. 66 in the PETSc Users Manual) 
    MatMult(J,unp1,b);
    VecAYPX(b, -1, un);

    // Solve for du
    KSPSolve(ksp, b, du);

    // Update u
    // TODO: Set unp1 = un + du
    //       and un = unp1
    // (Note: it turns out that could be done with just a single
    // vector un, but using two vectors makes the time increment
    // clearer)
    VecWAXPY(unp1, 1, un, du);
    VecCopy(unp1, un);

    // Store off middle value for output
    if (haveMiddle)
    {
      PetscScalar uMid;
      VecGetValues(unp1, 1, &Imid, &uMid);
      uMidValues[iTime+1] = uMid;
      timeValues[iTime+1] = (iTime+1)*dt;
    }

    // Output iteration info to screen
    int iters;
    KSPGetIterationNumber(ksp, &iters);
    double normDU;
    VecNorm(du, NORM_2, &normDU);
    PetscPrintf(PETSC_COMM_WORLD, "time = %7.4f, iters = %4d, |du| = %8.5g\n", (iTime+1)*dt, iters, normDU);
    
  }

  KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);

  // Stop timer and output timing
  double endTime = MPI_Wtime();
  PetscPrintf(PETSC_COMM_WORLD, "compute time = %10.4f\n", endTime - startTime);
  
  // Output middle temperature
  // TODO: Consider using command line options to set a filename instead of umiddle.dat (not required)
  if (haveMiddle)
  {
    writeValues(numTimeSteps+1, timeValues, uMidValues, name);
  }
  
  // Clean up
  VecDestroy(&un);
  VecDestroy(&unp1);
  VecDestroy(&du);
  VecDestroy(&b);
  MatDestroy(&J);
  KSPDestroy(&ksp);
  if (uMidValues) delete uMidValues;
  if (timeValues) delete timeValues;
  
  PetscFinalize();

}

// Write values to a file
void writeValues(int n, const double *t, const PetscScalar *u, const char fileName[])
{
  std::ofstream myFile;
  myFile.open(fileName);

  for (int i = 0; i < n; ++i)
  {
    char buf[256];
    char pattern[] = "%8.5f\t%9.6f";
    sprintf(buf, pattern, t[i], u[i]);
    myFile << buf << std::endl;
  }

  myFile.close();

}
  
  
