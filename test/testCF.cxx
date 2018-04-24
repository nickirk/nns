/*
 * testCF.cxx
 *
 *  Created on: Oct 26, 2017
 *      Author: guther
 */

#include <iostream>
#include "../src/FermionicHamiltonian.hpp"
#include "../src/CostFunction.hpp"
#include "../src/EnergyCF.hpp"
#include "../src/Determinant.hpp"

int main(){
	int numSites{3}, numStates{2*numSites};;
	double U{4}, t{-1};
	FermionicHamiltonian modelHam = generateFermiHubbard(numStates,U,t);
	EnergyCF eCF(modelHam);
	detType D1(numStates, false);
	detType D2 = D1;
	create(D1,0);
	create(D1,2);
	create(D2,0);
	create(D2,4);
	coeffType coeff = Eigen::VectorXd::Zero(2);
	coeff[0] = 1.0;
	std::vector<coeffType > TwoCoeffs(2,coeff);
	std::vector<detType > HFArr(2);
	HFArr[0] = D1;
	HFArr[1] = D2;
	std::vector<State> HF(1,State(HFArr,TwoCoeffs));
	std::vector<coeffType > dE = eCF.nabla(HF);
	std::cout << "Energy " << eCF.calc(HF) << std::endl;
	std::cout << "Derivative " << dE[0] << std::endl;
}

