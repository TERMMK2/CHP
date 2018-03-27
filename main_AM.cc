#include "fonction.h"
#include<iostream>
#include<cmath>
#include<vector>

int main(int argc, char * argv[])
{
  const int nx=10;
  const int ny=10;
  int N = nx*ny;
  int Me, Np, i1, iN;
  double beta,gamma,alpha,err,norm_r;
  //double alpha, beta, gamma, dx, dy;
  std::vector<std::vector<double> > A;
  A.resize(N);
  std::vector<double> b,d,p,r,r_1,x,w,D1,D2,D3,D4,D5,val_lock;
  b.resize(N); d.resize(N); p.resize(N); r.resize(N); r_1.resize(N); w.resize(N); x.resize(N); D1.resize(N); D2.resize(N); D3.resize(N); D4.resize(N); D5.resize(N);

  //Initialisation de beta, gamma, alpha et err
  beta=0. ; gamma=0.; alpha=0.; err=0.;

  //Initialisation des diagonales de la matrice A
  Diag_init(nx,ny,D1,D2,D3,D4,D5); 

  //Initialisation de x et de b
  for (int i=0 ; i<N; i++)
    {
      x[i] = 1.;
      b[i] = 1.;
    }

  //Initialisation de MPI
  MPI_Status status;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &Np); // get totalnodes
  MPI_Comm_rank(MPI_COMM_WORLD, &Me);

  charge(N,Np,Me,i1,iN);


  //Calcul de w
  prodMVC(argc,argv,D1,D2,D3,D4,D5,x,w,nx,ny);
  //Envoi de w
  if (Me==0)
      MPI_Bcast(&w[0],N,MPI_DOUBLE,Me,MPI_COMM_WORLD);

  //Calcul de r dans proc 0
  if (Me==0)
    {
      for (int i=0; i<N; i++)
	{
	  r[i] = b[i]-w[i];
	}
    }
  //Envoi de r
  if (Me==0)
    MPI_Bcast(&r[0],N,MPI_DOUBLE,Me,MPI_COMM_WORLD);
  
  //Calcul de norm_r
  norm_r=0.;
  //Envoi de norm_r
  if (Me==0)
    MPI_Bcast(&norm_r,1,MPI_DOUBLE,Me,MPI_COMM_WORLD);

  
  do while (norm_r > err)
       {
  	 //Calcul de d
  	 prodMVC(argc,argv,D1,D2,D3,D4,D5,x,w,nx,ny);
	 
  	 //Calcul de alpha
  	 alpha=0.;
	 
  	 //Calcul de x
  	 x = x+alpha*p;

  	 //Calcul de r_1
  	 r_1 = r-alpha*d;

  	 //Calcul de beta
  	 beta=0.;

  	 //Calcul de p
  	 p = r_1+gamma*p;
	 
  	 //Calcul de r_1
  	 r = r_1;  
       }


  if (Me==0)
    {
      for (int i=0; i<N; i++)
	{
	  std::cout << i << " " << x[i] << std::endl;
	}
    }

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
