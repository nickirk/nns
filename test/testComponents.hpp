/*
 * testComponents.hpp
 *
 *  Created on: Jul 27, 2018
 *      Author: guther
 */

#ifndef TEST_TESTCOMPONENTS_HPP_
#define TEST_TESTCOMPONENTS_HPP_

#include "../src/NNWLib.hpp"
// we also supply the default systems with the components
#include "defaultSystem.hpp"
// use ADAM to solve with EnergyEs
void solveADAM(networkVMC::Parametrization<> &para, networkVMC::Sampler &msampler, networkVMC::Hamiltonian const &H);
// use the full sampler with the direct parametrization
void testDeterministicFullSampling(networkVMC::SpinConfig const &sC, networkVMC::Hamiltonian const &H);
// check the adjacencyList of a lattice hamiltonian
void testAdj(networkVMC::LatticeHamiltonian const &test);


#endif /* TEST_TESTCOMPONENTS_HPP_ */
