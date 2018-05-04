#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <memory>
#include <vector>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "fonction.h"

#ifndef M_PI
   #define M_PI 3.141592653589793238462643383279
#endif



class Laplacian2D // pas fini de modifier
{

  protected: // Les attributs de la classe
    double _x_min, _x_max, _y_min, _y_max, _h_x, _h_y, _a, _deltaT;
    int _Nx, _Ny;
    std::vector<std::vector<double> > _LapMatloc; // matrice creuse du laplacien
    std::vector<double> _floc; // vecteur source _f qui prend les données de _sol(i) pour calculer _sol(i+1)
    std::vector<double> _solloc; // vecteur solution U local

    std::string _CL_bas, _CL_haut, _CL_gauche, _CL_droite;
    double _Val_CL_bas, _Val_CL_haut, _Val_CL_gauche, _Val_CL_droite;

    int _Me,_Np;

    std::string _Source;
    
    std::string _save_all_file;

    std::string _save_points_file;
    int _number_saved_points;
    std::vector<std::vector <double> > _saved_points;



  public: // Méthodes et opérateurs de la classe


    Laplacian2D();// Constructeur : Initialiser _x_min, _x_max, _y_min; _y_max; _N; _h; _LapMat; _x; _y et _sol.
    virtual ~Laplacian2D();


    void Initialize(double x_min, double x_max, double y_min, double y_max, int Nx, int Ny, double a, double deltaT, int Me, int Np, std::string Source, std::string save_all_file, std::string _save_points_file, int number_saved_points,std::vector<std::vector <double> > saved_points );

    void InitializeCL(std::string CL_bas, std::string CL_haut, std::string CL_gauche, std::string CL_droite, double Val_CL_bas, double Val_CL_haut, double Val_CL_gauche, double Val_CL_droite);

    void UpdateCL(int num_it);

    virtual void InitializeMatrix() = 0;

    void InitializeCI(double CI);

    virtual void IterativeSolver(int nb_iterations) = 0;   // Résout le système _LapMat * _sol = _f avec un solveur itératif.

    void SaveSol(std::string name_file); // Écrit les solutions dans le fichier "name_file".

    virtual void ConditionsLimites(int num_it) = 0;

  };

class EC_ClassiqueP : public Laplacian2D 
{
  //L'utilisation d'une unique classe fille vient de notre code utilisé dans notre TER sur lequel nous nous sommes basé pour écrire celui-ci. Même si elle n'est pas utile nous ne voulions pas la suprimer pour éviter de perdre du temps sur des choses peu utiles.


  public:
    void InitializeMatrix();
    void IterativeSolver(int nb_iterations);
    void ConditionsLimites(int num_it);
};
