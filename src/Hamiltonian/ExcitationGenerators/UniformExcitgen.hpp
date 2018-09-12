/*
 * UniformExcitgen.hpp
 *
 *  Created on: Jun 15, 2018
 *      Author: Lauretta Schwarz, guther
 */

#ifndef SRC_HAMILTONIAN_EXCITATIONGENERATORS_UNIFORMEXCITGEN_HPP_
#define SRC_HAMILTONIAN_EXCITATIONGENERATORS_UNIFORMEXCITGEN_HPP_

#include "ExcitationGenerator.hpp"
#include "ProbUpdater.hpp"

namespace networkVMC {

/**
 * \class UniformExcitgen
 * \brief Creates random excitations with a constant probability distribution
 *
 * Well suited for markov chain sampling, this ExcitationGenerator implementation picks a random
 * connected basis vector uniformly from a subset of possible excitations, with the subset being chosen
 * randomly with some dynamic probability distribution. These subsets are single excitations, same-spin
 * double excitations and opposite-spin double excitations.
 */
class UniformExcitgen: public ClonableExcitgen<UniformExcitgen> {
public:
	/// \param HF reference basis vector to initialize some internal parameters
	explicit UniformExcitgen(detType const &HF);
	virtual ~UniformExcitgen();
	virtual detType generateExcitation(detType const &source, double &pGen) ;
	virtual double getExcitationProb(detType const &source, detType const &target);

	// update the biases pDoubles/pParallel
	virtual void updateBiases();
private:
	// have separate methods for generating single/double excits
	detType genSingleExcitation(detType const &source, std::vector<int> const &source_orbs,
			std::vector<int> &holes, std::vector<int> &particles, double &pgen, double pdoubnew,
			std::vector<int> noccs, std::vector<int> nunoccs) const;

	detType genDoubleExcitation(detType const &source, std::vector<int> const &source_orbs,
			std::vector<int> &holes, std::vector<int> &particles, double &pgen, double pdoubnew,
			std::vector<int> noccs, std::vector<int> nunoccs) const;

	/// Maintains updating of pDoubles and pParallel
	ProbUpdater pBiasGen;
    /// probability of generating a double excitation
    double pDoubles;
    /// probability of generating a same-spin double excitation
    double pParallel;
    // which part of the excitation generator
    bool linExact = false;
    bool partExact = true;
};

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_UNIFORMEXCITGEN_HPP_ */
