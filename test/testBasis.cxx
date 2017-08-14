#include <stdio.h>
#include <iostream>
#include "../src/Determinant.hpp"
#include "../src/Basis.hpp"
using namespace std;
int main(){
  Basis basis(4,2);
  cout << "size= " << basis.getSize() << endl;
  for (int i=0; i < basis.getSize(); ++i){
    Determinant det;
    det = basis.getDetByIndex(i);
    vector<int> pos(det.getOccupiedPositions());
    for (int j=0; j < pos.size(); ++j){
      cout << pos[j] << "," ;
    }
    cout << endl;
  }
  return 0;
}
