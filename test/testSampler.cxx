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
