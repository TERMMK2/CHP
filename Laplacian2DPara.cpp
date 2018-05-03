#include "Laplacian2DPara.h"

using namespace std ;


//Constructeur :
Laplacian2D::Laplacian2D()
{}
//Destructeur :
Laplacian2D::~Laplacian2D()
{}

void Laplacian2D::Initialize(double x_min, double x_max, double y_min, double y_max, int Nx, int Ny, double a, double deltaT, int Me, int Np, string Source)
{
  // // On  initialise les constantes connues de tous les processeurs.

  _x_min = x_min;
  _y_min = y_min;
  _x_max = x_max;
  _y_max = y_max;
  _Nx = Nx;
  _Ny = Ny;
  _a = a;
  _deltaT = deltaT;
  _h_y = (y_min-y_max)/(Ny+1.);
  _h_x = (x_max-x_min)/(Nx+1.);
  _Me = Me;
  _Np = Np;
  _Source = Source;
}

void Laplacian2D::InitializeCI(double CI)
{
  // // On initialise le vecteur solution ici.

  int i1,iN;
  charge(_Nx*_Ny,_Np,_Me,i1,iN);

  _solloc.resize(iN-i1+1);
  for (int i=0; i<iN-i1+1;i++)
    {
      _solloc[i] = CI;
    }

  //On met _floc de la bonne taille ici  :
  _floc.resize(iN-i1+1);
}

void Laplacian2D::InitializeCL(std::string CL_bas, std::string CL_haut, std::string CL_gauche, std::string CL_droite, double Val_CL_bas, double Val_CL_haut, double Val_CL_gauche, double Val_CL_droite)
{
  // // On initialise les condition limites ici. La configuration quelque peu redondante avec Laplacian2D::Initialize vient d'une ancienne version du code de notre TER dans laquelle on initialisait toutes ces valeures dans le main et où on ne souhaitait pas avoir trop d'argument dans la méthode Initialize.

  _CL_bas = CL_bas;
  _CL_haut = CL_haut;
  _CL_gauche = CL_gauche;
  _CL_droite = CL_droite;
  _Val_CL_bas = Val_CL_bas;
  _Val_CL_haut = Val_CL_haut;
  _Val_CL_gauche = Val_CL_gauche;
  _Val_CL_droite = Val_CL_droite;
}

void Laplacian2D::UpdateCL(int num_it)
{
  //Cette méthode de classe nous permet de faire une condition au limites de Neumann variable au cours du temps. C'est avec cela que nous avons fait notre test pour vérifier notre code.

  double t = num_it*_deltaT;
  if (t <= 50.)
    {
      _Val_CL_gauche = 10000.*t;
    }
  else
    {
      _Val_CL_gauche = -9000.*(t-50.) + 500000.;
    }
}

void EC_ClassiqueP::InitializeMatrix()
{
  // // On initialise la matrice penta diagonale ici.



  vector<vector <double> > LapMat;
  LapMat.resize(5);
  int N = _Nx*_Ny;

  double alpha = 1 + 2*_a*_deltaT/(_h_x*_h_x) + 2*_a*_deltaT/(_h_y*_h_y);
  double beta = -_a*_deltaT/(_h_x*_h_x);
  double gamma = -_a*_deltaT/(_h_y*_h_y);

  // double alpha = 20;
  // double beta = 4;
  // double gamma = 5;

  for (int i=0; i<5; i++)
    {
      LapMat[i].resize(N);
    }



  for (int i=0; i<N; i++)
    {
      LapMat[0][i]=gamma;
      LapMat[1][i]=beta;
      LapMat[2][i]=alpha;
      LapMat[3][i]=beta;
      LapMat[4][i]=gamma;

      if (i%_Nx==0)
        LapMat[1][i]=0.;

      if (i%_Nx==_Nx-1) //hahaha les parenthèses mdr
        LapMat[3][i]=0.;

      if(i<_Nx)
        LapMat[0][i]=0.;

      if(i>(_Ny-1)*_Nx-1)
        LapMat[4][i]=0.;
    }

  if (_CL_gauche == "Neumann" or _CL_gauche == "Neumann_non_constant")
    {
      for (int i = 0 ; i < _Ny; i++)
        {
          LapMat[2][_Nx*i] += beta;  //Bord gauche
        }
    }

  if (_CL_droite == "Neumann")
    {
      for (int i = 0 ; i < _Ny; i++)
        {
          LapMat[2][_Nx*(i+1)-1] += beta; //Bord droit
        }
    }

  if (_CL_haut == "Neumann")
    {
      for (int i = 0; i < _Nx ; i++)
        {
          LapMat[2][i] += gamma; //Bord haut
        }
    }

  if (_CL_bas == "Neumann")
    {
      for (int i = 0; i < _Nx ; i++)
        {
          LapMat[2][(_Ny-1)*_Nx+i] += gamma; //Bord bas
        }
    }
  //--------------------------------------------------------------------
  //On découpe la matrice en vecteur local

  _LapMatloc.resize(5);
  int i1,iN;
  charge(_Nx*_Ny,_Np,_Me,i1,iN);

  for (int i=0; i<5; i++)
    {
      _LapMatloc[i].resize(iN-i1+1);
      _LapMatloc[i] = vectorsplit(LapMat[i]);
    }

  // printvect(_LapMatloc[3]);

}

void EC_ClassiqueP::IterativeSolver (int nb_iterations)
{
  // // Cette méthode est au coeur de la résolution du problème et elle nous permet de générer le résultat. <- PEUT MIEUX FAIRE

  //c'est le bordel pour  afficher les solutions en des points particuliers mais je vais changer ça quand j'aurais le temps, en attendant ça marche. <- Je le fais demain, j'ai pas eu le temps.

  system("rm -Rf EC_ClassiqueP");
  system("mkdir -p ./EC_ClassiqueP");

  MPI_Status status;

  int i1,iN;
  charge(_Nx*_Ny,_Np,_Me,i1,iN);
  int Nloc = iN-i1 +1;


  /*
  string _save_points_file ="sol_points_para";
  system(("rm -Rf "+_save_points_file).c_str());
  system(("mkdir -p ./"+_save_points_file).c_str());

  int _number_saved_points = 6;
  vector<vector <double> > _saved_points;
  _saved_points.resize(_number_saved_points);
  for (int i=0; i<6; i++)
    {
      _saved_points[i].resize(2);
      _saved_points[i][0] = 0.001*i;
      _saved_points[i][1] = 0.0025;
    }

 

  vector<double> _sol;


  vector< shared_ptr<ofstream> > mes_flux;

  if(_Me ==0)
    {
      for (int i=0; i<_number_saved_points; i++)
	{ 
	  shared_ptr<ofstream> flux(new ofstream);
	  flux->open(_save_points_file+"/point_"+to_string(i)+".txt", ios::out);
	  mes_flux.push_back(flux);
	}
    }
*/


  for( int i=0 ; i<=nb_iterations ; i++)
    {
      EC_ClassiqueP::SaveSol("EC_ClassiqueP/sol_it_"+to_string(i)+".vtk"); // -> a changer pour enregistrer la solution de chaque proc puis ensuite faire un truc pour reconbiner
      /*
      if(_Me ==0)
	{
	  
	  _sol.resize(_Nx*_Ny);
	  for (int i=0; i<=iN; i++)
	    {
	      _sol[i] = _solloc[i];
	    } 

	  for (int he=1; he<_Np ; he++)
	    {
	      int he_i1, he_iN;
	      charge(_Nx*_Ny,_Np,he,he_i1,he_iN);
	      vector<double> sol_temp;
	      sol_temp.resize(he_iN-he_i1+1);


	      MPI_Recv(&sol_temp[0],he_iN-he_i1+1,MPI_DOUBLE,he,100*he,MPI_COMM_WORLD, &status);

	      for (int i=he_i1; i<=he_iN ; i++)
		{
		  //cout<<_Me<<" "<<i<<" "<<sol_temp[i]<<endl;;
		  _sol[i] = sol_temp[i-he_i1];
		}
	    }

	  char* truc = new char;
	  for (int j=0; j<_number_saved_points; j++)
	    {
	      double truc_b1, truc_b2;
	      truc_b1 = i*_deltaT;
	      int pos = floor((_saved_points[j][0]/_h_x) + _Nx*floor(_saved_points[j][1]/_h_y));
	      truc_b2 = _sol[pos] ;
	      sprintf(truc, "%f  %f", truc_b1, truc_b2);
	      mes_flux[j]->write(truc,16);
	      mes_flux[j]->write("\n",1);
	    }
	}
      else
	{
	  MPI_Send(&_solloc[0],iN-i1+1,MPI_DOUBLE,0,100*_Me,MPI_COMM_WORLD);
	  }*/

      EC_ClassiqueP::ConditionsLimites(i);  // -> a changer pour faire des vecteurs locaux
      //--------------------------------------------------------------------------
      //Prise en compte du terme source :
      for (int k=0 ; k < Nloc ;k++)
	{
	  int num = i1 + k;
	  double x = (num%_Nx)*_h_x;
	  double y = (num/_Nx)*_h_y; 

	  if (_Source == "non")
	    {
	      _floc[k] = _solloc[k];
	    }
	  
	  if (_Source == "polynomial")
	    {
	       
	      _floc[k] = _solloc[k] + 2*(y - y*y + x -x*x);
	    }
	  if(_Source == "trigonometrique")
	    {
	      _floc[k] = _solloc[k] + sin(x) +cos(y) ;
	    }
	}
      
      //------------------------------------------------------------------------
      int kmax = _Nx*_Ny +100; //Pour une matrice de taille n le GC met max n étapes en théorie, comme on veut être sûr qu'il converge 
      
      _solloc = CGPara(_LapMatloc,_floc,_solloc,0.0001,kmax,_Nx,_Ny);

      if (_Me == 0)//Barre de chargement
	{
	  int i_barre;
	  int p = floor((((double)i)/((double)nb_iterations))*100);
	  printf( "[" );
	  for(i_barre=0;i_barre<=p;i_barre+=2) printf( "*" );
	  for (;i_barre<100; i_barre+=2 ) printf( "-" );
	  printf( "] %3d %%", p );
	  
	  for(i_barre=0;i_barre<59;++i_barre) printf( "%c", 8 );
	  
	  fflush(stdout );
	}


    }

  /*
  if(_Me ==0)
    {
      for (int i=0; i<_number_saved_points; i++)
	{ 
	  mes_flux[i]->close();
	}

    }
  */

  if (_Me == 0)//Barre de chargement
    { printf("\n");}


}

void Laplacian2D::SaveSol(string name_file)
{
  // // Cette méthode est celle qui nous permet d'enregistrer la solution sous forme de fichier lisible par paraview. Pour cela, elle envoie tout les vecteur locaux de la solution vers le processeur 0, qui va se charger de reformer le vecteur solution global puis d'écrire le résultat dans le bon format.

  MPI_Status status;

  int i1,iN;
  charge(_Nx*_Ny,_Np,_Me,i1,iN);
  if (_Me == 0)
    {
      vector<double> sol;
      sol.resize(_Nx*_Ny);
      for (int i=0; i<=iN; i++)
	{
	  sol[i] = _solloc[i];
	}

      for (int he=1; he<_Np ; he++)
	{
	  int he_i1, he_iN;
	  charge(_Nx*_Ny,_Np,he,he_i1,he_iN);
	  vector<double> sol_temp;
	  sol_temp.resize(he_iN-he_i1+1);


	  MPI_Recv(&sol_temp[0],he_iN-he_i1+1,MPI_DOUBLE,he,100*he,MPI_COMM_WORLD, &status);

	  for (int i=he_i1; i<=he_iN ; i++)
	    {
	      sol[i] = sol_temp[i-he_i1];
	    }
	}

      using namespace std;
      ofstream mon_flux;
      mon_flux.open(name_file, ios::out);
      mon_flux << "# vtk DataFile Version 3.0" << endl
	       << "cell" << endl
	       << "ASCII" << endl
	       << "DATASET STRUCTURED_POINTS" << endl
	       << "DIMENSIONS " << _Nx << " " << _Ny << " 1" << endl
	       << "ORIGIN 0 0 0" << endl
	       << "SPACING " + to_string((_x_max-_x_min)/_Nx)+ " " + to_string((_y_max-_y_min)/_Ny) +" 1" << endl
	       << "POINT_DATA " << _Nx*_Ny << endl
	       << "SCALARS sample_scalars double" << endl
	       << "LOOKUP_TABLE default" << endl;

      for(int i=_Ny-1; i>=0; i--)
	{
	  for(int j=0; j<_Nx; j++)
	    {
	      mon_flux << sol[j + i*_Nx] << " ";
	    }
	  mon_flux << endl;
	}

      mon_flux.close();
    }
  else
    {
      MPI_Send(&_solloc[0],iN-i1+1,MPI_DOUBLE,0,100*_Me,MPI_COMM_WORLD);
    }

}


void EC_ClassiqueP::ConditionsLimites(int num_it)
{
  // // Cette méthode nous permet de mettre à jours le terme source à chaque itération pour prendre en compte les effets des conditions limites.

  int i1,iN;
  double gamma = -_a*_deltaT/(_h_y*_h_y);
  double beta = -_a*_deltaT/(_h_x*_h_x);
  charge(_Nx*_Ny,_Np,_Me,i1,iN);



  if (_CL_haut == "Dirichlet") //Condition de température en haut
    {
      for (int j = 0; j < _Nx ; j++)
        {
	  if ((j<=iN) and (j>=i1))
	    _solloc[j-i1] = _solloc[j-i1]-gamma*_Val_CL_haut;
        }
    }
  if (_CL_bas == "Dirichlet") //Condition de température en bas
    {
      for (int j = 0; j < _Nx ; j++)
        {
	  if ((_Nx*(_Ny -1)+ j<=iN) and (_Nx*(_Ny -1)+ j>=i1))
	    _solloc[_Nx*(_Ny -1)+ j -i1] = _solloc[_Nx*(_Ny -1)+ j -i1]-gamma*_Val_CL_bas;
        }
    }

  if (_CL_gauche == "Dirichlet")  //Condition de température à gauche
    {
      for (int i = 0; i < _Ny; i++)
        {
	  if ((i*_Nx<=iN) and (i*_Nx>=i1))
	    _solloc[i*_Nx -i1] = _solloc[i*_Nx -i1]-beta*_Val_CL_gauche;
        }
    }

  if (_CL_droite == "Dirichlet") //Condition de température à droite
    {
      for (int i = 0; i < _Ny; i++)
        {
	  if (((i+1)*_Nx - 1 <= iN) and ((i+1)*_Nx - 1 >= i1))
	    _solloc[(i+1)*_Nx - 1 -i1] = _solloc[(i+1)*_Nx - 1 -i1]-beta*_Val_CL_droite;
        }
    }

  if (_CL_haut == "Neumann") //Condition de flux en haut
    {
      for (int j = 0; j < _Nx ; j++)
        {
	  if ((j<=iN) and (j>=i1))
	    _solloc[j -i1] = _solloc[j -i1]-gamma*_Val_CL_haut*_h_y;
        }
    }
  if (_CL_bas == "Neumann") //Condition de flux en bas
    {
      for (int j = 0; j < _Nx ; j++)
        {
	  if ((_Nx*(_Ny -1)+ j<=iN) and (_Nx*(_Ny -1)+ j>=i1))
	    _solloc[_Nx*(_Ny -1)+ j -i1] = _solloc[_Nx*(_Ny -1)+ j -i1]-gamma*_Val_CL_bas*_h_y;
        }
    }

  if (_CL_gauche == "Neumann")  //Condition de flux à gauche
    {
      for (int i = 0; i < _Ny; i++)
        {
	  if ((i*_Nx<=iN) and (i*_Nx>=i1))
	    _solloc[i*_Nx -i1] = _solloc[i*_Nx -i1]-beta*_Val_CL_gauche*_h_x;
        }
    }

  if (_CL_gauche == "Neumann_non_constant")  //Condition de flux à gauche
    {
      Laplacian2D::UpdateCL(num_it);
      for (int i = 0; i < _Ny; i++)
        {
	  if ((i*_Nx<=iN) and (i*_Nx>=i1))
	    _solloc[i*_Nx -i1] = _solloc[i*_Nx -i1]-beta*_Val_CL_gauche*_h_x;
        }
    }

  if (_CL_droite == "Neumann") //Condition de flux à droite
    {
      for (int i = 0; i < _Ny; i++)
        {
	  if (((i+1)*_Nx - 1 <= iN) and ((i+1)*_Nx - 1 >= i1))
	    _solloc[(i+1)*_Nx - 1 -i1] = _solloc[(i+1)*_Nx - 1 -i1]-beta*_Val_CL_droite*_h_x;
        }
     }


  // if (_Source == "trigonometrique")
  //   {
  //     for (int j = 0; j < _Nx ; j++)
  //       {
  // 	  double x
  // 	  if ((j<=iN) and (j>=i1))
  // 	    _solloc[j-i1] = _solloc[j-i1]-gamma*_Val_CL_haut;
  //       }




  //   }

}
