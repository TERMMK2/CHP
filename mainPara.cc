
#include "Laplacian2DPara.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

using namespace std;

int main()
{
  int Nx = 300;
  int Ny = 2;
  double xmin = 0.;
  double xmax = 0.01;
  double ymin = 0.;
  double ymax = 0.005;
  double a = 1./(1500*1000);
  double deltaT = 0.05;
  double tfinal = 100.;
  std::vector<double> CI;
  string CL_bas = "Neumann"; // "Neumann" , "Dirichlet"
  string CL_haut = "Neumann";
  string CL_gauche = "Neumann_non_constant"; //Peut valoir "Neumann_non_constant"
  string CL_droite = "Neumann";
  double Val_CL_bas = 0; //Flux si CL_bas == "Neumann", Température si CL_bas == "Dirichlet"
  double Val_CL_haut = 0;
  double Val_CL_gauche = 0; //Mettre 0 si CL_gauche == "Neumann_non_constant"
  double Val_CL_droite = 0;
  int nb_iterations = int(ceil(tfinal/deltaT));
  string Equation = "EC_ClassiqueP";

  CI.resize(Nx*Ny);
  for(int j=0; j < Ny; j++)
  {
    for(int i=0; i < Nx; i++)
    {
      CI[i + j*Nx] = 293.;
    }
  }

  Laplacian2D *Lap;

  Lap = new EC_ClassiqueP();
  Lap->Initialize(xmin,xmax,ymin,ymax,Nx,Ny,a,deltaT);
  Lap->InitializeCI(CI);
  Lap->InitializeCL(CL_bas, CL_haut, CL_gauche, CL_droite, Val_CL_bas, Val_CL_haut, Val_CL_gauche, Val_CL_droite);
  Lap->InitializeMatrix();
  Lap->IterativeSolver(nb_iterations);

  return 0;
}
