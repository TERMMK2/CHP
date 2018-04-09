#include "fonction.h"


void charge(int N, int Np, int Me, int& i1, int& iN)
{
  int q = N/Np;
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

// int prodMV(int argc, char * argv[], std::vector<std::vector<double> > A, std::vector<double> x)
// {
//
//   int N;
//   N = x.size();
//   int i1, iN, j1, jN, r;
//   int Me, Np;
//   std::vector<double> y,yloc;
//   y.resize(N); yloc.resize(N);
//
//
//   for (int i=0; i<N; i++)
//     {
//       y[i] = 0;
//     }
//
//   MPI_Status status;
//   MPI_Init(&argc,&argv);
//   MPI_Comm_size(MPI_COMM_WORLD, &Np); // get totalnodes
//   MPI_Comm_rank(MPI_COMM_WORLD, &Me);
//
//   charge(N,Np,Me,i1,iN);
//
//   std::cout << "Je suis le proc " << Me << " je commence à : " << i1 << " et je termine à : " << iN << std::endl;
//
//
//   //Calcul local de Y
//
//   yloc.resize(iN-i1+1);
//
//   for (int i=i1; i<iN+1; i++)
//     {
//       yloc[i-i1] = 0;
//       for (int j=0; j<N; j++)
// 	{
// 	  yloc[i-i1] = yloc[i-i1] + A[i][j]*x[j];
// 	}
//     }
//
//   if(Me!=0){
//     MPI_Send(&yloc[0],iN-i1+1,MPI_DOUBLE,0,100,MPI_COMM_WORLD);
//   }
//   else
//     {
//       for (int i=i1; i<iN+1; i++){
// 	y[i] = yloc[i-i1];
//       }
//       for (int he=1; he<Np; he++){
//
// 	charge(N,Np,he,j1,jN);
//
// 	yloc.resize(jN-j1+1);
//
// 	MPI_Recv(&yloc[0],jN-j1+1,MPI_DOUBLE,he,100,MPI_COMM_WORLD, &status);
//
// 	for (int i=j1; i<jN+1; i++){
// 	  y[i] = yloc[i-j1];
// 	}
//       }
//     }
//
//
//   if (Me==0){
//     std::cout<<"Valeur du produit AX ="<<std::endl;
//     for (int i=0; i<N; i++)
//       {
// 	std::cout << i << " " << y[i] << std::endl;
//       }
//   }
//
//   MPI_Finalize();
//   return 0;
// }
//
// //-------------------------------------------------------------
//

std::vector<double> prodMVC(std::vector<std::vector<double> > Aloc, std::vector<double> xloc, int nx, int ny)
{
  int N, i1, iN, Me, Np;
  N = nx*ny;
  std::vector<double> x_haut, x_bas, yloc; //x_haut contient les éléments voisins du dessus
  x_bas.resize(nx);
  x_haut.resize(nx);

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &Np); // get totalnodes
  MPI_Comm_rank(MPI_COMM_WORLD, &Me);

  charge(N,Np,Me,i1,iN);

  int Nloc = iN-i1+1;

  if (Me%2 == 0)
  {
    if (Me != Np-1)
    {
      MPI_Send(&xloc[Nloc-nx],nx,MPI_DOUBLE,Me+1,100,MPI_COMM_WORLD);
    }
    if (Me != 0)
    {
      MPI_Send(&xloc[0],nx,MPI_DOUBLE,Me-1,100,MPI_COMM_WORLD);
    }
  }
  else
  {
    MPI_Recv(&x_haut[0],nx,MPI_DOUBLE,Me-1,100,MPI_COMM_WORLD, &status);
    if (Me != Np-1)
    {
      MPI_Recv(&x_bas[0],nx,MPI_DOUBLE,Me+1,100,MPI_COMM_WORLD, &status);
    }
  }

  if (Me%2 == 1)
  {
    if (Me != Np-1)
    {
      MPI_Send(&xloc[Nloc-nx],nx,MPI_DOUBLE,Me+1,100,MPI_COMM_WORLD);
    }
    MPI_Send(&xloc[0],nx,MPI_DOUBLE,Me-1,100,MPI_COMM_WORLD);
  }
  else
  {
    if (Me != 0)
    {
      MPI_Recv(&x_haut[0],nx,MPI_DOUBLE,Me-1,100,MPI_COMM_WORLD, &status);
    }
    if (Me != Np-1)
    {
      MPI_Recv(&x_bas[0],nx,MPI_DOUBLE,Me+1,100,MPI_COMM_WORLD, &status);
    }
  }

  //Calcul local de Y

  yloc.resize(iN-i1+1);

  for (int i =0 ; i < iN-i1+1; i++)
  {
    yloc[i] = 0;
    yloc[i] += Aloc[2][i]*xloc[i];
    if (Me == 0) //Pas de voisin en haut sur la première ligne
    {
      if (i == 0) //pas de voisin à gauche
      {
        yloc[i] += Aloc[3][i]*xloc[i+1];
      }
      else if (i == iN - i1) //voisin à droite dans x_bas
      {
        yloc[i] += Aloc[1][i]*xloc[i-1] + Aloc[3][i]*x_bas[0];
      }
      else //voisins gauche et droite dans xloc
      {
        yloc[i] += Aloc[1][i]*xloc[i-1] + Aloc[3][i]*xloc[i+1];
      }


      if (i < iN-i1+1-nx) //voisin du bas dans xloc
      {
        yloc[i] += Aloc[4][i]*xloc[i+nx];
      }
      else //voisin du bas dans x_bas
      {
        yloc[i] += Aloc[4][i]*x_bas[i - (iN-i1-nx+1)];
      }


      if (i >= nx) //voisin du haut dans xloc à partir de la deuxième ligne
      {
        yloc[i] += Aloc[0][i]*xloc[i-nx];
      }
    }


    else if (Me == Np-1) //Pas de voisin en bas sur la dernière ligne
    {
      if (i == iN - i1) //pas de voisin à droite
      {
        yloc[i] += Aloc[1][i]*xloc[i-1];
      }
      else if (i == 0) //voisin à gauche dans x_haut
      {
        yloc[i] += Aloc[3][i]*xloc[i+1] + Aloc[1][i]*x_haut[nx-1];
      }
      else //voisins gauche et droite dans xloc
      {
        yloc[i] += Aloc[1][i]*xloc[i-1] + Aloc[3][i]*xloc[i+1];
      }


      if (i < iN-i1+1-nx) //voisin du bas dans xloc
      {
        yloc[i] += Aloc[4][i]*xloc[i+nx];
      }


      if (i >= nx) //voisin du haut dans xloc à partir de la deuxième ligne
      {
        yloc[i] += Aloc[0][i]*xloc[i-nx];
      }
      else //voisin du haut dans x_haut
      {
        yloc[i] += Aloc[0][i]*x_haut[i];
      }
    }


    else //Cas général
    {
      if (i == 0) //voisin à gauche dans x_haut
      {
        yloc[i] += Aloc[1][i]*x_haut[nx-1] + Aloc[3][i]*xloc[i+1];
      }
      else if (i == iN - i1) //voisin à droite dans x_bas
      {
        yloc[i] += Aloc[1][i]*xloc[i-1] + Aloc[3][i]*x_bas[0];
      }
      else //voisins gauche et droite dans xloc
      {
        yloc[i] += Aloc[1][i]*xloc[i-1] + Aloc[3][i]*xloc[i+1];
      }


      if (i < iN-i1+1-nx) //voisin du bas dans xloc
      {
        yloc[i] += Aloc[4][i]*xloc[i+nx];
      }
      else //voisin du bas dans x_bas
      {
        yloc[i] += Aloc[4][i]*x_bas[i - (iN-i1-nx+1)];
      }


      if (i >= nx) //voisin du haut dans xloc
      {
        yloc[i] += Aloc[0][i]*xloc[i-nx];
      }
      else //voisin du haut dans x_haut
      {
        yloc[i] += Aloc[0][i]*x_haut[i];
      }
    }


  }
  return yloc;
}

double dot(std::vector<double> uloc, std::vector<double> vloc)
{
  int Nloc = uloc.size();
  double y,yloc;

  //Calcul d'une partie de y
  yloc = 0;
  for (int i=0; i<Nloc; i++)
  {
    yloc += uloc[i]*vloc[i];
  }

  MPI_Allreduce(&yloc,&y,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  return y;
}

std::vector<double> vectorsplit(std::vector<double> u)
{
  int N, Np, Me, i1, iN;
  std::vector<double> uloc;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &Np); // get totalnodes
  MPI_Comm_rank(MPI_COMM_WORLD, &Me);
  N = u.size();
  charge(N,Np,Me,i1,iN);
  uloc.resize(iN - i1 + 1);
  for (int i = 0; i <= iN - i1; i++)
  {
    uloc[i] = u[i + i1];
  }
  return uloc;
}

void Diag_init(int nx, int ny, std::vector<std::vector<double> >& A)
{

  double alpha, beta, gamma, dx, dy;
  int N = nx*ny;
  dx = 1./(nx+1);
  dy = 1./(ny+1);
  alpha = 2./(dx*dx)+2./(dy*dy);
  beta = -1./(dx*dx);
  gamma = -1./(dy*dy);
  A.resize(5);
  for(int i = 0; i < 5; i++)
  {
    A[i].resize(N);
  }

  for (int i=0; i<N; i++)
    {
      A[0][i]=5;
      A[1][i]=4;
      A[2][i]=20;
      A[3][i]=4;
      A[4][i]=5;

  //     if (i%nx==0)
	// A[1][i]=0.;
  //
  //     if (i%(nx-1)==0)
	// A[3][i]=0.;
  //
  //     if(i<nx)
	// A[0][i]=0.;
  //
  //     if(i>(ny-1)*nx-1)
	// A[4][i]=0.;
    }
}

std::vector<double> CGPara (std::vector<std::vector<double> > Aloc, std::vector<double> bloc, std::vector<double> x0loc , double err, int kmax, int nx, int ny)
{

  int i1, iN, Nloc, k, Me, Np;

  double norm_r, nr_carre, nr2_carre;
  std::vector<double> wloc, rloc, r_1loc, ploc, dloc, xloc;

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &Np); // get totalnodes
  MPI_Comm_rank(MPI_COMM_WORLD, &Me);

  charge(nx*ny,Np,Me,i1,iN);

  Nloc = iN-i1 +1;

  wloc.resize(Nloc); rloc.resize(Nloc); ploc.resize(Nloc); dloc.resize(Nloc);

  xloc = x0loc;

  wloc = prodMVC(Aloc,xloc,nx,ny);

  for (int i = 0; i < Nloc; i++)
  {
    rloc[i] = bloc[i] - wloc[i];
  }
  ploc = rloc;

  k = 0;

  nr_carre = dot(rloc,rloc);
  norm_r = sqrt(nr_carre);


  while ((norm_r > err) and (k<kmax))
    {
      dloc = prodMVC(Aloc,ploc,nx,ny);

      double alpha = nr_carre/dot(ploc,dloc);

      for (int i = 0; i < Nloc; i++)
      {
        xloc[i] += alpha*ploc[i];
        rloc[i] -= alpha*dloc[i]; //rk devient r(k+1)
      }

      nr2_carre = dot(rloc,rloc); //norme de r(k+1)

      double beta = nr2_carre/nr_carre;

      for (int i = 0; i < Nloc; i++)
      {
        ploc[i] = rloc[i] + beta*ploc[i];
      }

      nr_carre = nr2_carre;
      norm_r = sqrt(nr_carre);

      k++;
    }


  return xloc;

}

void printvect(std::vector<double> uloc)
{
  int Me, Np;
  MPI_Comm_rank(MPI_COMM_WORLD, &Me);
  MPI_Comm_size(MPI_COMM_WORLD, &Np); // get totalnodes
  for (int i =0; i < Np; i++)
    {
      if (Me == i)
      {
        std::cout << "Me = " << i << std::endl;
        for (int j = 0; j < uloc.size(); j++)
        {
          std::cout << uloc[j] << std::endl;
        }
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  if (Me == Np-1)
  {
    std::cout << "fin affichage" << std::endl;
  }
}
