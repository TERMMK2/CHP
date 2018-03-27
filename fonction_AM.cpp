#include "fonction.h"


void charge(int N, int Np, int Me, int& i1, int& iN)
{
  int q = floor((double)N/(double)Np);
  int r = N - q*Np;

  if (Me < r)
    {
      i1 = (q+1)*Me;
      iN = (q+1)*(Me+1)-1;
    }
  else
    {
      i1 = q*Me + r;
      iN = q*(Me+1) + r - 1;
    }
}

//---------------------------------------------------------------

int prodMV(int argc, char * argv[], std::vector<std::vector<double> > A, std::vector<double> x)
{
  
  int N;
  N = x.size();
  int i1, iN, j1, jN, r;
  int Me, Np;
  std::vector<double> y,yloc;
  y.resize(N); yloc.resize(N);

  
  for (int i=0; i<N; i++)
    {
      y[i] = 0;
    }

  MPI_Status status;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &Np); // get totalnodes
  MPI_Comm_rank(MPI_COMM_WORLD, &Me);
  
  charge(N,Np,Me,i1,iN);

  std::cout << "Je suis le proc " << Me << " je commence à : " << i1 << " et je termine à : " << iN << std::endl;


  //Calcul local de Y
  
  yloc.resize(iN-i1+1);
  
  for (int i=i1; i<iN+1; i++)
    {
      yloc[i-i1] = 0;
      for (int j=0; j<N; j++)
	{
	  yloc[i-i1] = yloc[i-i1] + A[i][j]*x[j];
	}
    }
  
  if(Me!=0){
    MPI_Send(&yloc[0],iN-i1+1,MPI_DOUBLE,0,100,MPI_COMM_WORLD);
  }
  else
    {
      for (int i=i1; i<iN+1; i++){
	y[i] = yloc[i-i1];
      }
      for (int he=1; he<Np; he++){

	charge(N,Np,he,j1,jN);

	yloc.resize(jN-j1+1);

	MPI_Recv(&yloc[0],jN-j1+1,MPI_DOUBLE,he,100,MPI_COMM_WORLD, &status);

	for (int i=j1; i<jN+1; i++){
	  y[i] = yloc[i-j1];
	}
      }
    }

  
  if (Me==0){
    std::cout<<"Valeur du produit AX ="<<std::endl;
    for (int i=0; i<N; i++)
      {
	std::cout << i << " " << y[i] << std::endl;
      }
  }

  MPI_Finalize(); 
  return 0;
}

//-------------------------------------------------------------

int prodMVC(int argc, char * argv[], std::vector<double> D1, std::vector<double> D2, std::vector<double> D3, std::vector<double> D4, std::vector<double> D5, std::vector<double> x, std::vector<double>& y, int nx, int ny)
{
  
  int N;
  N = x.size();
  int i1, iN, j1, jN, r;
  int Me, Np;
  std::vector<double> yloc;
  yloc.resize(N);

  
  for (int i=0; i<N; i++)
    {
      y[i] = 0;
    }

  MPI_Status status;
  //MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &Np); // get totalnodes
  MPI_Comm_rank(MPI_COMM_WORLD, &Me);
  
  charge(N,Np,Me,i1,iN);

  //std::cout << "Je suis le proc " << Me << " je commence à : " << i1 << " et je termine à : " << iN << std::endl;


  //Calcul local de Y
  
  yloc.resize(iN-i1+1);
  
  for (int i=i1; i<iN+1; i++){ 

    if (i<nx){ //La diagonale D1 n'intervient pas dans le calcul
  
      if (i==0)
	yloc[i-i1] = (D3[i]*x[i] + D4[i]*x[i+1]) + D5[i]*x[i+nx];

      else
      yloc[i-i1] = (D2[i]*x[i-1] + D3[i]*x[i] + D4[i]*x[i+1]) + D5[i]*x[i+nx];
    }
    
    else if (i>((ny-1)*nx-1)){ //La diagonale D5 n'intervient pas dans le calcul

      if (i==N-1)
	yloc[i-i1] = D1[i]*x[i-nx] + (D2[i]*x[i-1] + D3[i]*x[i]);

      else
      yloc[i-i1] = D1[i]*x[i-nx] + (D2[i]*x[i-1] + D3[i]*x[i] + D4[i]*x[i+1]);
    }

    else
      yloc[i-i1] = D1[i]*x[i-nx] + (D2[i]*x[i-1] + D3[i]*x[i] + D4[i]*x[i+1]) + D5[i]*x[i+nx];
      
  }    
   
  
  if(Me!=0){
    MPI_Send(&yloc[0],iN-i1+1,MPI_DOUBLE,0,100,MPI_COMM_WORLD);
  }
  else
    {
      for (int i=i1; i<iN+1; i++){
	y[i] = yloc[i-i1];
      }
      for (int he=1; he<Np; he++){

	charge(N,Np,he,j1,jN);

	yloc.resize(jN-j1+1);

	MPI_Recv(&yloc[0],jN-j1+1,MPI_DOUBLE,he,100,MPI_COMM_WORLD, &status);

	for (int i=j1; i<jN+1; i++){
	  y[i] = yloc[i-j1];
	}
      }
    }

  
  // if (Me==0){
  //   std::cout<<"Valeur du produit AX ="<<std::endl;
  //   for (int i=0; i<N; i++)
  //     {
  // 	std::cout << i << " " << y[i] << std::endl;
  //     }
  // }

  //MPI_Finalize(); 
  return 0;
}


//-----------------------------------------------------------------

void Diag_init(int nx, int ny, std::vector<double>& D1, std::vector<double>& D2, std::vector<double>& D3, std::vector<double>& D4, std::vector<double>& D5)
{

  double alpha, beta, gamma, dx, dy;
  int N = nx*ny;
  
  dx = 1./(nx+1);
  dy = 1./(ny+1);
  alpha = 2./(dx*dx)+2./(dy*dy);
  beta = -1./(dx*dx);
  gamma = -1./(dy*dy);

  for (int i=0; i<N; i++)
    {
      D1[i]=gamma;
      D2[i]=beta;
      D3[i]=alpha; 
      D4[i]=beta;
      D5[i]=gamma;
      
      if (i%nx==0)
	D2[i]=0.;

      if (i%(nx-1)==0)
	D4[i]=0.;

      if(i<nx)
	D1[i]=0.;

      if(i>(ny-1)*nx-1)
	D5[i]=0.;
    }
}
