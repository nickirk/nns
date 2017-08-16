#include <stdio.h>
#include <iostream>
#include <Eigen/Eigenvalues>
#include "../src/Determinant.hpp"
#include "../src/Basis.hpp"
#include "../src/Hamiltonian.hpp"
using namespace std;
int main(){
  Basis basis(5,2);
  Hamiltonian modelHam(0.5, 0.2, 0.2, basis);
  cout << "Basis size= " << basis.getSize() << endl;
  cout << "Hamiltonian size= " << modelHam.getSize() << endl;
  cout << "print out Ham element" << endl;
  for (int i(0); i<modelHam.getSize(); ++i){
    for (int j(0); j< modelHam.getSize(); ++j){
      cout << modelHam(i,j) << ", " ;
    }
    cout << endl;
  }
  Eigen::MatrixXd H(Eigen::MatrixXd::Zero(modelHam.getSize(),modelHam.getSize())); 
  for (int i=0; i<modelHam.getSparseSize(); ++i){
    int row, col;
    double value(0.);
    modelHam.sparseAccess(i, row, col, value);
    H(row,col) = value;
  }
  cout << "Diag= " << H.eigenvalues()  << endl;
}

