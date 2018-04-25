#include <stdio.h>
#include <iostream>
#include <Eigen/Eigenvalues>
#include "../src/Hamiltonian/FermionicHamiltonian.hpp"
#include "../src/HilbertSpace/Basis.hpp"
#include "../src/HilbertSpace/Determinant.hpp"
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
  cout << "Hamiltonian size= " << modelHam.getNumOrbs() << endl;
  cout << "print out Ham element" << endl;
  for (int i(0); i<modelHam.getNumOrbs(); ++i){
    for (int j(0); j< modelHam.getNumOrbs(); ++j){
      cout << modelHam(i,j) << ", " ;
    }
    cout << endl;
  }
  for (int i(0); i<modelHam.getNumOrbs(); ++i){
    for (int j(0); j< modelHam.getNumOrbs(); ++j){
      cout << modelHam(basis.getDetByIndex(i), basis.getDetByIndex(j)) << ",";
    }
    cout << endl;
  }
}

