#include <vector>
#include <iostream>
#include <math.h>
#include "../src/Nnw.hpp"
#include "../src/Basis.hpp"
#include "../src/Determinant.hpp"
#include "../src/Hamiltonian.hpp"
using namespace std;
int main(){
  int numStates(7);
  int numEle(3);
  Basis basis(numStates,numEle);
  Hamiltonian modelHam(0.5, 0.2, 0.3, basis);
  cout << "Basis size= " << basis.getSize() << endl;
  cout << "Hamiltonian size= " << modelHam.getSize() << endl;
  cout << "print out Ham element" << endl;
  Eigen::MatrixXd H(Eigen::MatrixXd::Zero(basis.getSize(),basis.getSize()));
  for (int i(0); i<modelHam.getSize(); ++i){
    for (int j(0); j< modelHam.getSize(); ++j){
      H(i,j) = modelHam(i,j); 
      //cout << modelHam(i,j) << ", " ;
    }
  }
  cout << "eigenV= \n" << H.eigenvalues() << endl;
  vector<int> size_NNW = {numStates, 40,1};
  vector<detType> list;
  for (int i=0; i<basis.getSize(); ++i){
    list.push_back(basis.getDetByIndex(i));
    vector<int> pos=getOccupiedPositions(basis.getDetByIndex(i));
    for (int j=0; j<pos.size(); j++){
      cout << pos[j] << ",";
    }
    cout << endl;
  }
  NeuralNetwork NNW(size_NNW, modelHam, basis);
  while (true){
    NNW.train(list, 0.4); 
    cout << "energy= " << NNW.getEnergy() << endl;
  }
}
