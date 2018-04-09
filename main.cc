#include "fonction.h"
#include<iostream>
#include<cmath>
#include<vector>

int main(int argc, char * argv[])
{
  MPI_Status status;
  MPI_Init(&argc,&argv);
  int Me, Np;
  double err = 1.e-4;
  int kmax = 100000;

  const int nx=2;
  const int ny=3;
  int N = nx*ny;
  //double alpha, beta, gamma, dx, dy;
  std::vector<std::vector<double> > A, Aloc;
  std::vector<double> x,x0,bloc,xloc,x0loc;
  x.resize(N);
  x0.resize(N);


  for (int i=0 ; i<N; i++)
  {
    x[i] = 1.;
    x0[i] = 0.;
  }

  xloc = vectorsplit(x);
  x0loc = vectorsplit(x0);

  Diag_init(nx,ny,A); //Initialisation des diagonales de la matrice A
  Aloc.resize(5);
  for (int i = 0; i < 5; i++)
  {
    Aloc[i] = vectorsplit(A[i]);
  }
  bloc = prodMVC(Aloc,xloc,nx,ny); //Produit Matrice-Vecteur creux

  // if (Me == 0)
  // {
  //   std::cout << "b =" << std::endl;
  // }
  // printvect(bloc);

  xloc = CGPara (Aloc, bloc, x0loc, err, kmax, nx, ny);

  MPI_Comm_rank(MPI_COMM_WORLD, &Me);

  if (Me == 0)
  {
    std::cout << "x =" << std::endl;
  }
  printvect(xloc);

  MPI_Finalize();
  return 0;
}










 // for(int i=0; i<N; i++)
  //   {
  //     A[i].resize(N);
  //     x[i] = 1.;
  //     if (i<nx){ 
  // 	if (i==0){
  // 	  A[i][0] = D3[i];
  // 	  A[i][1] = D4[i];
  // 	  A[i][nx] = D5[i];
  // 	}
  // 	else{
  // 	  A[i][i-1] = D2[i];
  // 	  A[i][i] = D3[i];
  // 	  A[i][i+1] = D4[i];
  // 	  A[i][i+nx] = D5[i];
  // 	}
  //     }

  //     else if (i>((ny-1)*nx-1)){ 
  // 	if (i==N-1){
  // 	  A[i][N-2] = D2[i];
  // 	  A[i][N-1] = D3[i];
  // 	  A[i][N-1-nx] = D1[i];
  // 	}
  // 	else{
  // 	  A[i][i-1] = D2[i];
  // 	  A[i][i] = D3[i];
  // 	  A[i][i+1] = D4[i];
  // 	  A[i][i-nx] = D1[i];
  // 	}
  //     }

  //     else{
  // 	A[i][i-1] = D2[i];
  // 	A[i][i] = D3[i];
  // 	A[i][i+1] = D4[i];
  // 	A[i][i-nx] = D1[i];
  // 	A[i][i+nx] = D5[i];
  //     }

  //   }
  
  //prodMV(argc, argv, A, x); //Produit Matrice-Vecteur
