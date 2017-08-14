#include <stdio.h>
#include <iostream>
#include "../src/Determinant.hpp"
#include "../src/Basis.hpp"
#include "../src/Hamiltonian.hpp"
using namespace std;
int main(){
  Basis basis(10,4);
  Hamiltonian modelHam(1.5, 0.2, 0.2, basis);
  cout << "Basis size= " << basis.getSize() << endl;
  cout << "Hamiltonian size= " << modelHam.getSize() << endl;
}
