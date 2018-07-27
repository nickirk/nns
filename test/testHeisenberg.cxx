/*
 * testHeisenberg.cxx
 *
 *  Created on: Feb 5, 2018
 *      Author: guther
 */

#include "../src/Hamiltonian/HeisenbergHamiltonian.hpp"
#include <iostream>

#include "../src/HilbertSpace/Determinant.hpp"

using namespace networkVMC;

int main(){
	double J{2};
	HeisenbergHamiltonian test(J,4,4);
	std::cout << test.size() << std::endl;
	for(int i = 0; i < test.size(); ++i){
		std::cout << "Adjacent to " << i << ": ";
		for(auto j:test.adj(i)){
			std::cout << j << " ";
		}
		std::cout<<std::endl;
	}
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


