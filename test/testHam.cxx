#include <stdio.h>
#include <iostream>
#include <Eigen/Eigenvalues>
#include "../src/Hamiltonian/FermiHubbardHamiltonian.hpp"
#include "../src/HilbertSpace/FermionBasis.hpp"
#include "../src/HilbertSpace/Determinant.hpp"
#include "defaultSystem.hpp"
using namespace std;
using namespace networkVMC;
int main(){
  int numSites(8);
  int numStates(2*numSites);
  int numEle(8);
  int spinUp(4);
  int spinDown(4);
  SpinConfig spinConfig{spinUp, spinDown, numStates};
  //cout << "input number of hidden neurons=";
  //cin >> numHidden;
  bool readFromFile{false};
  double trainRate(1.5);
  FermionBasis basis(spinConfig);
  auto modelHam = generateDefaultHubbard(numSites);
  cout << "Basis size= " << basis.getSize() << endl;
  cout << "Hamiltonian size= " << modelHam.getNumOrbs() << endl;
  cout << "print out Ham element" << endl;
  for (int i(0); i<modelHam.getNumOrbs(); ++i){
    for (int j(0); j< modelHam.getNumOrbs(); ++j){
      cout << modelHam(basis.getDetByIndex(i), basis.getDetByIndex(j)) << ",";
    }
    cout << endl;
  }
}

