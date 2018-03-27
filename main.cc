#include "fonction.h"
#include<iostream>
#include<cmath>
#include<vector>

int main(int argc, char * argv[])
{
  const int nx=10;
  const int ny=10;
  int N = nx*ny;
  //double alpha, beta, gamma, dx, dy;
  std::vector<std::vector<double> > A;
  A.resize(N);
  std::vector<double> x,D1,D2,D3,D4,D5;
  x.resize(N); D1.resize(N); D2.resize(N); D3.resize(N); D4.resize(N); D5.resize(N);

 

  for (int i=0 ; i<N; i++)
    x[i] = 1.;

  Diag_init(nx,ny,D1,D2,D3,D4,D5); //Initialisation des diagonales de la matrice A
  prodMVC(argc,argv,D1,D2,D3,D4,D5,x,nx,ny); //Produit Matrice-Vecteur creux

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
