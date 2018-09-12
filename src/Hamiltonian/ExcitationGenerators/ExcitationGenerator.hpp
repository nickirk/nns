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

/**
 * \class ExcitationGenerator
 * \brief Abstract base class for excitation generation
 *
 * Interface for ExcitationGenerators, objects that can randomly generate basis vectors connected to a given one.
 * In addition, they can output the probability of generating a given basis vector from another given one.
 */
class ExcitationGenerator {
public:
	ExcitationGenerator(){};
	virtual ~ExcitationGenerator(){};

	/**
	 * \brief Virtual constructor (implemented with CRTP via ClonableExcitgen)
	 * \return pointer to a new ExcitationGenerator which is a copy of *this (with the same derived type)
	 */
	virtual ExcitationGenerator* clone()const =0;

	// we need to be able to
	//	a) generate a random excitation and

	/**
	 * \brief Create a random excitation of a given basis vector
	 * \param[in] source basis vector from which to excite
	 * \param[out] pGen probability to pick this excitation
	 * \return random coupled basis vector (excitation)
	 */
	virtual detType generateExcitation(detType const &source, double &pGen) = 0;
	// 	b) get the probability for an excitation from source to target

	/**
	 * \brief Get the probability to create a given excitation
	 * \param[in] source initial basis vector of the excitation generation
	 * \param[in] target final basis vector of the excitation generation
	 * \return probability to excite from source to target with generateExcitation()
	 */
	virtual double getExcitationProb(detType const &source, detType const &target) = 0;
	/// maintains the pDoubles/pParallel biases for our non-trivial excitgens
	class ProbUpdater;
	/** additional feature: update internal biases of excitation generation, this can make a huge difference
	* since this is not done on the same scope as generating excitations, we need to call this externally
	*/
	virtual void updateBiases(){};
};


/**
 * \class ClonableExcitgen
 * \brief CRTP class for cloning ExcitationGenerators
 * \tparam T derived type that inherits from this instance
 */

template<typename T>
class ClonableExcitgen: public ExcitationGenerator{
public:
	/// inherit the constructor of ExcitationGenerator
	using ExcitationGenerator::ExcitationGenerator;
	virtual ~ClonableExcitgen(){};

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
