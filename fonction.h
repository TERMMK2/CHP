#include<iostream>
#include<mpi.h>
#include<cmath>
#include<vector>

void charge (int N, int Np, int Me, int& i1, int& iN);

// int prodMV(int argc, char * argv[], std::vector<std::vector<double> > A, std::vector<double> x);

std::vector<double> prodMVC(std::vector<std::vector<double> > Aloc, std::vector<double> xloc, int nx, int ny);

double dot(std::vector<double> u, std::vector<double> v);

std::vector<double> vectorsplit(std::vector<double> u);

void Diag_init(int nx, int ny, std::vector<std::vector<double> >& A);

std::vector<double> CGPara (std::vector<std::vector<double> > Aloc, std::vector<double> bloc, std::vector<double> x0loc , double err, int kmax, int nx, int ny);

void printvect( std::vector<double> uloc);
