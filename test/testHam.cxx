#include <stdio.h>
#include <iostream>
#include <fstream>
#include <Eigen/Eigenvalues>
#include "../src/Hamiltonian/FermiHubbardHamiltonian.hpp"
#include "../src/HilbertSpace/FermionBasis.hpp"
#include "../src/HilbertSpace/Determinant.hpp"
#include "testComponents.hpp"
using namespace std;
using namespace networkVMC;

int main(){
  int numSites = 4;
  auto basis = generateDefaultBasis(numSites);
  auto modelHam = generateDefaultHubbard(numSites);
  testAdj(modelHam);
  cout << "Basis size= " << basis.size() << endl;
  cout << "Hamiltonian size= " << modelHam.size() << endl;
  cout << "print out Ham element" << endl;
  ofstream ofs;
  ofs.open("ham");
  for (int i(0); i<basis.size(); ++i){
    for (int j(0); j< basis.size(); ++j){
      cout << modelHam(basis.getDetByIndex(i), basis.getDetByIndex(j)) << ",";
      ofs << modelHam(basis.getDetByIndex(i), basis.getDetByIndex(j)) << ",";
    }
    cout << endl;
    ofs << endl;
  }
  ofs.close();
}

