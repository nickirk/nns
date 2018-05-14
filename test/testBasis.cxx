#include <stdio.h>
#include <iostream>

#include "../src/HilbertSpace/Basis.hpp"
#include "../src/HilbertSpace/Determinant.hpp"
//Basis.cxx and Basis.hpp have been tested and passed
using namespace std;
using namespace networkVMC;

int main(){
  int numSites(8);
  int numStates(2*numSites);
  int numEle(8);
  int spinUp(4);
  int spinDown(4);
  vector<int> spinConfig{spinUp, spinDown, numStates};
  Basis basis(spinConfig);
  cout << "size= " << basis.getSize() << endl;
  for (int i=0; i < basis.getSize(); ++i){
    detType det;
    det = basis.getDetByIndex(i);
    int index(basis.getIndexByDet(det));
    vector<int> pos(getOccupiedPositions(det));
    cout << "index= " << index << endl;
    cout << "i = " << i  << endl;
    for (int j=0; j < pos.size(); ++j){
      cout << pos[j] << "," ;
    }
    cout << endl;
  }
  return 0;
}
