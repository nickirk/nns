#include <stdio.h>
#include <iostream>
#include "../src/Determinant.hpp"
#include "../src/Basis.hpp"
//Basis.cxx and Basis.hpp have been tested and passed
using namespace std;

int main(){
  Basis basis(8,4);
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
