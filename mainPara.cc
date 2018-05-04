
#include "Laplacian2DPara.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <chrono>

using namespace std;

// COMMENTER UN PEU !
// GENRE METTRE DES TITRES.


int main(int argc, char * argv[])
{
  MPI_Status status;
  MPI_Init(&argc,&argv);
  int Me, Np;
  MPI_Comm_size(MPI_COMM_WORLD, &Np); // get totalnodes
  MPI_Comm_rank(MPI_COMM_WORLD, &Me);

  int Nx = 300;
  int Ny = 300;

  double xmin = 0.;
  double xmax = 1.;
  double ymin = 0.;
  double ymax = 1.;

  double a = 1.;
  double deltaT = 0.005;
  double tfinal = 1.;

  string CL_bas = "Dirichlet"; // "Neumann" , "Dirichlet"
  string CL_haut = "Dirichlet";
  string CL_gauche = "Dirichlet"; //Peut valoir "Neumann_non_constant"
  string CL_droite = "Dirichlet";

  double Val_CL_bas = 0; //Flux si CL_bas == "Neumann", Température si CL_bas == "Dirichlet"
  double Val_CL_haut = 0;
  double Val_CL_gauche = 0; 
  double Val_CL_droite = 0;




  int nb_iterations = int(ceil(tfinal/deltaT));
  string Equation = "EC_ClassiqueP";

  string Source = "polynomial"; //Peut prendre non, polynomial ou trigonometrique. 
  //Choisir trigonométrique met à jour les conditions limites directement.

  double CI = 0.1;

  
  string save_all_file = "Polyomial";
  
  string save_points_file = "non";
  int number_saved_points=6;
  vector<vector <double> > saved_points;
  saved_points.resize(number_saved_points);
  for (int i=0; i<6; i++)
    {
      saved_points[i].resize(2);
      saved_points[i][0] = 0.001*i;
      saved_points[i][1] = 0.0025;
    }

  
  Laplacian2D *Lap;
  Lap = new EC_ClassiqueP();


  Lap->Initialize(xmin,xmax,ymin,ymax,Nx,Ny,a,deltaT, Me, Np, Source, save_all_file, save_points_file, number_saved_points, saved_points);
  Lap->InitializeCI(CI);
  Lap->InitializeCL(CL_bas, CL_haut, CL_gauche, CL_droite, Val_CL_bas, Val_CL_haut, Val_CL_gauche, Val_CL_droite);
  Lap->InitializeMatrix();


  auto start = chrono::high_resolution_clock::now();
  Lap->IterativeSolver(nb_iterations);
  auto finish = chrono::high_resolution_clock::now();

  double t = chrono::duration_cast<chrono::microseconds>(finish-start).count();

  if(Me ==0)
    {
      cout << "Le prog a mis " << t*0.000001 << " secondes a s'effectuer" << endl;
    }
  MPI_Finalize();

  return 0;
}
