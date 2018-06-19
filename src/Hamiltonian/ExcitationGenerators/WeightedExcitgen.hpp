/*
 * WeightedExcitgen.hpp
 *
 *  Created on: Jun 15, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_EXCITATIONGENERATORS_WEIGHTEDEXCITGEN_HPP_
#define SRC_HAMILTONIAN_EXCITATIONGENERATORS_WEIGHTEDEXCITGEN_HPP_

//this excitation generator uses the matrix elements of a Hamiltonian to get
//generation probabilities

#include "ExcitationGenerator.hpp"

class Hamiltonian;

namespace networkVMC {

class WeightedExcitgen: public ExcitationGenerator {
public:
	WeightedExcitgen();
	virtual ~WeightedExcitgen();

	// interface for excitation generation
	virtual detType generateExcitation(detType const &source, double &pGen) const ;
	virtual double getExcitationProb(detType const &source, detType const &target) const;

private:
    // generate a single excitation
    detType generateSingleExcit();
    // generate a double excitation
    detType generateDoubleExcit();
    std::vector<int> WeightedExcitgen::pickBiasedElecs(std::vector<int> &elecs, detType const &source);
    double pParallel, pDoubles, pGen;

    // indicate whether all variables are filled
    bool bfilled;
    // number of spin orbitals
    int norbs;
    // number of electrons
    int nel;
    // number of alpha spin electrons
    int nalphaels;
    // number of beta spin electrons
    int nbetaels;
    // number of alpha spin holes
    int nalphaholes;
    // number of beta spin holes
    int nbetaholes;
    // number of pairs of alpha-alpha spin electrons
    int aa_elec_pairs;
    // number of pairs of beta-beta spin electrons
    int bb_elec_pairs;
    // number of pairs of alpha-beta spin electrons
    int ab_elec_pairs;

	Hamiltonian const *H;
};

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_WEIGHTEDEXCITGEN_HPP_ */
