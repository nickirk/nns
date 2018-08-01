/*
 * testHeisenberg.cxx
 *
 *  Created on: Feb 5, 2018
 *      Author: guther
 */

#include "../src/Hamiltonian/HeisenbergHamiltonian.hpp"
#include <iostream>
#include "testComponents.hpp"

#include "../src/HilbertSpace/Determinant.hpp"

using namespace networkVMC;

int main(){
	double J{-1.0};
	HeisenbergHamiltonian test(J,false,4,4);
	// test the Heisenberg Basis
	SpinConfig sC(2,2,4);
	Basis basis(sC,test);
	testBasis(basis);
	// test the adjacency of the Hamiltonian
	testAdj(test);
	// and then chech some matrix elements
	detType a(test.size(),false);
	std::cout << "FM energy: " << test(a,a) << std::endl;
	detType b = a;
	detType c =a;
	create(b,1);
	create(c,13);
	std::cout << "Coupling energy: " << test(b,c) << std::endl;
	create(c,6);
	std::cout << "Different spins: " << test(b,c) << std::endl;
}


