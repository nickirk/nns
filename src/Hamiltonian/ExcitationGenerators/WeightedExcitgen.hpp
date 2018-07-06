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
#include "ProbUpdater.hpp"
#include "ExcitmatType.hpp"
#include "../Hamiltonian.hpp"

namespace networkVMC {

class WeightedExcitgen: public clonableExcitgen<WeightedExcitgen> {
public:
	// constructed from a Hamiltonian to get the weights and a
	// determinant to use as a starting guess for the bias pDoubles
	WeightedExcitgen(Hamiltonian const &H_, detType const &HF);
	virtual ~WeightedExcitgen();

	// interface for excitation generation
	virtual detType generateExcitation(detType const &source, double &pGen);
	virtual double getExcitationProb(detType const &source, detType const &target);

	// update the biases pDoubles/pParallel
	virtual void updateBiases();
private:
    // generate a single excitation
    detType generateSingleExcit(detType const &source, double &hel);
    // can throw a NoExcitFound exception
    // generate a double excitation
    detType generateDoubleExcit(detType const &source, double &hel);
    // can throw a NoExcitFound exception

    // pick two electrons
    std::vector<int> pickBiasedElecs(
    		std::vector<int> &elecs, detType const &source);

    // compute the probability to generate an excitation moving the electrons
    // src to the orbitals tgt
    double calcPgen(detType const &source,std::vector<int> const &src,
    		std::vector<int> const &tgt);
    // calculate the generation probability for a single excitation
    double pGenSingleExcitCs(detType const &source, int srcOrb, int tgt);
    // calculate the generation probability for selecting holes in a doule excitation
    double pGenSelectDoubleHole(
    		detType const &source, std::vector<int> const &srcOrbs,
			int const &orb_pair, double &cum_sum, int const &tgt);

    // Initialize all the internal variables
    void constructClassCount(detType const &source);

    ProbUpdater pBiasGen;
    double pParallel, pDoubles, pgen;
    // auxiliary matrix for internal excitation communication
    ExcitmatType excitmat;
    // orbitals occupied in a source
    std::vector<int> sourceOrbs;
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
