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
// use a trial WF parametrized RBM with metropolis sampling
void testTrialMetropolis(networkVMC::SpinConfig const &sC, networkVMC::Hamiltonian const &H);
// check the adjacencyList of a lattice hamiltonian
void testAdj(networkVMC::LatticeHamiltonian const &test);
// check the basis generation
void testBasis(networkVMC::Basis const &basis);


// use ADAM to solve
template<typename F, typename coeffType>
void solveEs(networkVMC::Parametrization<F, coeffType> &para, networkVMC::Sampler<coeffType> &msampler,
		networkVMC::Hamiltonian const &H, networkVMC::Solver<F, coeffType> &solver){
	networkVMC::EnergyEs<F, coeffType> eCF(H,-1);
	networkVMC::Trainer<F, coeffType> etr(para,msampler,solver,eCF,H);
	networkVMC::detType HF;
	for(int i = 1; i < 4; ++i){
		etr.train();
		if(i==1) HF = etr.getState().det(0);
		std::cout << "On iteration " << i << std::endl;
		std::cout << "Current energy: " << etr.getE()<<std::endl;
		std::cout << "Current state: " << std::endl;
		auto pars = para.pars();
		std::cout << pars(0) << std::endl;
		std::cout << para.getCoeff(HF) << std::endl;
	}
}

template<typename F, typename coeffType>
void solveADAM(networkVMC::Parametrization<F, coeffType> &para, networkVMC::Sampler<coeffType> &msampler,
		networkVMC::Hamiltonian const &H){
	networkVMC::ADAM<F, coeffType> solver(0.01);
    solveEs(para,msampler,H,solver);
}

template<typename F, typename coeffType>
void solveSRec(networkVMC::Parametrization<F, coeffType> &para, networkVMC::Sampler<coeffType> &msampler, networkVMC::Hamiltonian const &H){
	networkVMC::StochasticReconfiguration<F, coeffType> solver(para,0.01);
	solveEs(para,msampler,H,solver);
}



#endif /* TEST_TESTCOMPONENTS_HPP_ */
