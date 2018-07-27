#include <stdio.h>
#include <iostream>

#include "../src/HilbertSpace/Basis.hpp"
#include "../src/HilbertSpace/Determinant.hpp"
//Basis.cxx and Basis.hpp have been tested and passed
using namespace std;
using namespace networkVMC;

int main(){
  int numSites(4);
  int numStates(2*numSites);
  int numEle(numSites);
  int spinUp(numSites/2);
  int spinDown(numSites/2);
  SpinConfig spinConfig{spinUp, spinDown, numStates};
  Basis basis(spinConfig);
  cout << "size= " << basis.size() << endl;
  for (int i=0; i < basis.size(); ++i){
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
