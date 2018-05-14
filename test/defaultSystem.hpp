/*
 * defaultSystem.hpp
 *
 *  Created on: May 4, 2018
 *      Author: guther
 */

#ifndef TEST_DEFAULTSYSTEM_HPP_
#define TEST_DEFAULTSYSTEM_HPP_

#include "../src/Hamiltonian/FermionicHamiltonian.hpp"
#include "../src/HilbertSpace/Basis.hpp"

networkVMC::FermionicHamiltonian generateDefaultHubbard(int numSites);
networkVMC::Basis generateDefaultBasis(int numSites);


#endif /* TEST_DEFAULTSYSTEM_HPP_ */
