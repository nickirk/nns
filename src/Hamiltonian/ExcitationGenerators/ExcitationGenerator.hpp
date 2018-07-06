/*
 * ExcitationGenerator.hpp
 *
 *  Created on: Jun 15, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_EXCITATIONGENERATORS_EXCITATIONGENERATOR_HPP_
#define SRC_HAMILTONIAN_EXCITATIONGENERATORS_EXCITATIONGENERATOR_HPP_

#include "../../utilities/TypeDefine.hpp"
#include <memory>

namespace networkVMC {

// Base class for excitation generation
class ExcitationGenerator {
public:
	ExcitationGenerator(){};
	virtual ~ExcitationGenerator(){};

	// we want polymorphic copy, it is implemented with CRTP (see below)
	virtual ExcitationGenerator* clone()const =0;

	// we need to be able to
	//	a) generate a random excitation and
	virtual detType generateExcitation(detType const &source, double &pGen) = 0;
	// 	b) get the probability for an excitation from source to target
	virtual double getExcitationProb(detType const &source, detType const &target) = 0;
	// maintains the pDoubles/pParallel for our non-trivial excitgens
	class ProbUpdater;
	// additional feature: update internal biases of excitation generation, this can make a huge difference
	// since this is not done on the same scope as generating excitations, we need to call this externally
	virtual void updateBiases(){};
};


// CRTP class for cloning

template<typename T>
class ClonableExcitgen: public ExcitationGenerator{
public:
	// inherit the constructor
	using ExcitationGenerator::ExcitationGenerator;
	virtual ~ClonableExcitgen(){};

	// clone functionality for ExcitationGenerators
	// only usable in direct child classes
	virtual ExcitationGenerator* clone() const{
			return new T{static_cast<T const&>(*this)};
	}
};

// utilities for general selection of electrons/orbitals

int countForbiddenOrbs(std::vector<int> const &spins, std::vector<int> nunoccs);
std::vector<int> pickElecPair(std::vector<int> const &source_orbs, std::vector<int> &spin, int &elecpairs);

int pickOrbA(detType const &source, std::vector<int> const &spins, std::vector<int> noccs,
		std::vector<int> nunoccs, int norbs, int nel, int nforbiddenorbs, int &nexcita,
		int &aspin, bool &baorbfail);

int pickOrbB(detType const &source, std::vector<int> const &spins, std::vector<int> noccs,
		std::vector<int> nunoccs, int norbs, int nel, int nforbiddenorbs, int aorb,
		int aspin, int &nexcitb, int &bspin, int &nexcitotherway);

} /* namespace networkVMC */



#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_EXCITATIONGENERATOR_HPP_ */
