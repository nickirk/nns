/*
 * defaultSystem.hpp
 *
 *  Created on: May 4, 2018
 *      Author: guther
 */

#ifndef TEST_DEFAULTSYSTEM_HPP_
#define TEST_DEFAULTSYSTEM_HPP_

#include "../src/Hamiltonian/FermiHubbardHamiltonian.hpp"
#include "../src/HilbertSpace/FermionBasis.hpp"

networkVMC::FermiHubbardHamiltonian generateDefaultHubbard(int numSites);
networkVMC::FermionBasis generateDefaultBasis(int numSites);

#endif /* TEST_DEFAULTSYSTEM_HPP_ */
