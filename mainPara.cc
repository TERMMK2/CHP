
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

  int Nx = 500;
  int Ny = 10;

  double xmin = 0.;
  double xmax = 0.04;
  double ymin = 0.;
  double ymax = 0.005;

  double a = 1./(1500.*1000.); //Mettre 1. si on fait les cas tests de l'énoncé et 1./(1500.*1000.) si on veut comparer avec notre TER.
  double deltaT = 0.05;
  double tfinal = 100;

  string CL_bas = "Neumann"; // "Neumann" , "Dirichlet"
  string CL_haut = "Neumann";
  string CL_gauche = "Neumann_non_constant"; //Peut valoir "Neumann_non_constant"
  string CL_droite = "Neumann";

  double Val_CL_bas = 0; //Flux si CL_bas == "Neumann", Température si CL_bas == "Dirichlet"
  double Val_CL_haut = 0;
  double Val_CL_gauche = 0; 
  double Val_CL_droite = 0;




  int nb_iterations = int(ceil(tfinal/deltaT));

  string Source = "non"; //Peut prendre non, polynomial, trigonometriqueou instationnaire. 
  //Choisir trigonométrique met à jour les conditions limites automatiquement.

  double CI = 293.;

  
  string save_all_file = "TER"; //Mettre non si on ne souhaite pas enregistrer la solution globale au cours du temps sous une forme lisible par paraview
  
  string save_points_file = "points_TER"; //Mettre non si l'on 
  int number_saved_points=3;
  vector<vector <double> > saved_points;
  saved_points.resize(number_saved_points);
  for (int i; i<number_saved_points; i++)
    {
      saved_points[i].resize(2);
      saved_points[i][0] = 0.002*i;
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
