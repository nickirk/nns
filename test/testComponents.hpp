/*
 * testComponents.hpp
 *
 *  Created on: Jul 27, 2018
 *      Author: guther
 */

#ifndef TEST_TESTCOMPONENTS_HPP_
#define TEST_TESTCOMPONENTS_HPP_
#include <iostream>
#include "../src/NNWLib.hpp"
// we also supply the default systems with the components
#include "defaultSystem.hpp"
// use the full sampler with the direct parametrization
void testDeterministicFullSampling(networkVMC::SpinConfig const &sC, networkVMC::Hamiltonian const &H);
// use the RBM with metropolis sampling
void testRBMMetropolis(networkVMC::SpinConfig const &sC, networkVMC::Hamiltonian const &H);
// check the adjacencyList of a lattice hamiltonian
void testAdj(networkVMC::LatticeHamiltonian const &test);
// check the basis generation
void testBasis(networkVMC::Basis const &basis);


// use ADAM to solve
template<typename T>
void solveEs(networkVMC::Parametrization<T> &para, networkVMC::Sampler &msampler,
		networkVMC::Hamiltonian const &H, networkVMC::Solver<T> &solver){
	networkVMC::EnergyEs eCF(H,-1);
	networkVMC::Trainer<T> etr(para,msampler,solver,eCF,H);
	for(int i = 1; i < 4000; ++i){
		etr.train();
		std::cout << "On iteration " << i << std::endl;
		std::cout << "Current energy: " << etr.getE()<<std::endl;
	}
}

template<typename T>
void solveADAM(networkVMC::Parametrization<T> &para, networkVMC::Sampler &msampler,
		networkVMC::Hamiltonian const &H){
	networkVMC::ADAM<T> solver(0.01);
    solveEs(para,msampler,H,solver);
}

template<typename T>
void solveSRec(networkVMC::Parametrization<T> &para, networkVMC::Sampler &msampler, networkVMC::Hamiltonian const &H){
	networkVMC::StochasticReconfiguration<T> solver(para,0.01);
	solveEs(para,msampler,H,solver);
}



#endif /* TEST_TESTCOMPONENTS_HPP_ */
