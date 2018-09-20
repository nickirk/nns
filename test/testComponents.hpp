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


// use ES to solve (then uses ADAM)
void solveEs(networkVMC::Parametrization &para, networkVMC::Sampler &msampler,
		networkVMC::Hamiltonian const &H, networkVMC::Solver &solver);

void solveADAM(networkVMC::Parametrization &para, networkVMC::Sampler &msampler,
		networkVMC::Hamiltonian const &H);

void solveSRec(networkVMC::Parametrization &para, networkVMC::Sampler &msampler, networkVMC::Hamiltonian const &H);



#endif /* TEST_TESTCOMPONENTS_HPP_ */
