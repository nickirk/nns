/*
 * ExcitationGenerator.hpp
 *
 *  Created on: Jun 15, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_EXCITATIONGENERATORS_EXCITATIONGENERATOR_HPP_
#define SRC_HAMILTONIAN_EXCITATIONGENERATORS_EXCITATIONGENERATOR_HPP_

#include "../../utilities/TypeDefine.hpp"

namespace networkVMC {

// Base class for excitation generation
class ExcitationGenerator {
public:
	ExcitationGenerator();
	virtual ~ExcitationGenerator(){};

	// we need to be able to
	//	a) generate a random excitation and
	virtual detType generateExcitation(detType const &source, double &pGen) const = 0;
	// 	b) get the probability for an excitation from source to target
	virtual double getExcitationProb(detType const &source, detType const &target) const = 0;
};

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
