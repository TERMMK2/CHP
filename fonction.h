#include<iostream>
#include<mpi.h>
#include<cmath>
#include<vector>

void charge (int N, int Np, int Me, int& i1, int& iN);

int prodMV(int argc, char * argv[], std::vector<std::vector<double> > A, std::vector<double> x);

int prodMVC(int argc, char * argv[], std::vector<double> D1, std::vector<double> D2, std::vector<double> D3, std::vector<double> D4, std::vector<double> D5, std::vector<double> x, int nx, int ny);

void Diag_init(int nx, int ny, std::vector<double>& D1, std::vector<double>& D2, std::vector<double>& D3, std::vector<double>& D4, std::vector<double>& D5);

