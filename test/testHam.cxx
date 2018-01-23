#include <stdio.h>
#include <iostream>
#include <Eigen/Eigenvalues>
#include "../src/Determinant.hpp"
#include "../src/Basis.hpp"
#include "../src/Hamiltonian.hpp"
using namespace std;
int main(){
  int numSites(8);
  int numStates(2*numSites);
  int numEle(8);
  int spinUp(4);
  int spinDown(4);
  vector<int> spinConfig{spinUp, spinDown, numStates};
  //cout << "input number of hidden neurons=";
  //cin >> numHidden;
  bool readFromFile{false};
  double trainRate(1.5);
  cout << "input training rate=";
  cin >> trainRate;
  Basis basis(spinConfig);
  FermionicHamiltonian modelHam(numStates);
  cout << "Basis size= " << basis.getSize() << endl;
  cout << "Hamiltonian size= " << modelHam.getSize() << endl;
  cout << "print out Ham element" << endl;
  for (int i(0); i<modelHam.getSize(); ++i){
    for (int j(0); j< modelHam.getSize(); ++j){
      cout << modelHam(i,j) << ", " ;
    }
    cout << endl;
  }
  for (int i(0); i<modelHam.getSize(); ++i){
    for (int j(0); j< modelHam.getSize(); ++j){
      cout << modelHam(basis.getDetByIndex(i), basis.getDetByIndex(j)) << ",";
    }
    cout << endl;
  }
  /*
  cout << modelHam.lowerPos(0)<<'\t'<<modelHam.upperPos(0)<<'\n';
  int a,b;
  double c;
  for(int i(0);i<modelHam.getSparseSize();++i){
    modelHam.sparseAccess(i,a,b,c);
    cout << a << '\t'<<b<<'\t'<<c<<'\n';
  }
  */
  Eigen::MatrixXd H(Eigen::MatrixXd::Zero(modelHam.getSize(),modelHam.getSize())); 
  for (int i=0; i<modelHam.getSparseSize(); ++i){
    int row, col;
    double value(0.);
    modelHam.sparseAccess(i, row, col, value);
    H(row,col) = value;
  }
  cout << "Diag= " << H.eigenvalues()  << endl;
}

