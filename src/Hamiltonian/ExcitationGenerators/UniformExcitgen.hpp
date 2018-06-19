/*
 * UniformExcitgen.hpp
 *
 *  Created on: Jun 15, 2018
 *      Author: Lauretta Schwarz, guther
 */

#ifndef SRC_HAMILTONIAN_EXCITATIONGENERATORS_UNIFORMEXCITGEN_HPP_
#define SRC_HAMILTONIAN_EXCITATIONGENERATORS_UNIFORMEXCITGEN_HPP_

#include "ExcitationGenerator.hpp"

namespace networkVMC {

class UniformExcitgen: public ExcitationGenerator {
public:
	explicit UniformExcitgen(detType const &HF);
	virtual ~UniformExcitgen();
	virtual detType generateExcitation(detType const &source, double &pGen) const ;
	virtual double getExcitationProb(detType const &source, detType const &target) const;
private:
	// have separate methods for generating single/double excits
	detType genSingleExcitation(detType const &source, std::vector<int> const &source_orbs,
			std::vector<int> &holes, std::vector<int> &particles, double &pgen, double pdoubnew,
			std::vector<int> noccs, std::vector<int> nunoccs) const;

	detType UniformExcitgen::genDoubleExcitation(detType const &source, std::vector<int> const &source_orbs,
			std::vector<int> &holes, std::vector<int> &particles, double &pgen, double pdoubnew,
			std::vector<int> noccs, std::vector<int> nunoccs) const;

    // probability of generating a double excitation
    double pDoubles;
    // probability of generating a parallel spin double excitation
    double pParallel;
    // dynamic adjustment of the probabilities
    bool bbiasSd;
    bool bbiasPo;
    // which part of the excitation generator
    bool linExact = false;
    bool partExact = true;
    // set values for the generation probabilities based on the number
    // of electrons and holes
    void setProbabilities(detType example_det);
    // initialise these parameters
    void setProbabilitiesBias(detType example_det, int exflag);
    // for updating these parameters
    void checkProbabilities(excitStore const &excitation, double hel, int nspawns);
    //void check_probabilities(double hel, int nspawns);
    void adjustProbabilities();
};

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_UNIFORMEXCITGEN_HPP_ */
