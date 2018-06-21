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

class UniformExcitgen: public clonableExcitgen<UniformExcitgen> {
public:
	explicit UniformExcitgen(detType const &HF);
	virtual ~UniformExcitgen();
	virtual detType generateExcitation(detType const &source, double &pGen) ;
	virtual double getExcitationProb(detType const &source, detType const &target);
private:
	// have separate methods for generating single/double excits
	detType genSingleExcitation(detType const &source, std::vector<int> const &source_orbs,
			std::vector<int> &holes, std::vector<int> &particles, double &pgen, double pdoubnew,
			std::vector<int> noccs, std::vector<int> nunoccs) const;

	detType genDoubleExcitation(detType const &source, std::vector<int> const &source_orbs,
			std::vector<int> &holes, std::vector<int> &particles, double &pgen, double pdoubnew,
			std::vector<int> noccs, std::vector<int> nunoccs) const;

	ProbUpdater pBiasGen;
    // probability of generating a double excitation
    double pDoubles;
    // probability of generating a parallel spin double excitation
    double pParallel;
    // which part of the excitation generator
    bool linExact = false;
    bool partExact = true;
};

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_UNIFORMEXCITGEN_HPP_ */
