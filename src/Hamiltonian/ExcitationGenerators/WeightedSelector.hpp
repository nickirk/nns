/*
 * WeightedSelector.hpp
 *
 *  Created on: Jun 19, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_EXCITATIONGENERATORS_WEIGHTEDSELECTOR_HPP_
#define SRC_HAMILTONIAN_EXCITATIONGENERATORS_WEIGHTEDSELECTOR_HPP_

#include "../../utilities/TypeDefine.hpp"
// auxiliary excitation generation functionality:
// given a Hamiltonian and a determinant, we can select electrons based
// on the Hamiltonian matrix elements
// it is essentially a wrapper to create cumulative sums

namespace networkVMC{

class Hamiltonian;

class WeightedSelector {
public:
	// this is created on the fly to select some parts/holes from a determinant
	// called source, weighted with a Hamiltonian H
	WeightedSelector(Hamiltonian const &H_, detType const &source_):
		H(H_),source(source_),norbs(source_.size()){};
	virtual ~WeightedSelector();
    // select a hole for a single excitation (hel is the matrix element of that excitation)
    int selectSingleHole(int src, double &pgen, double &hel);
    // select two holes based on a cumulative list
    int selectDoubleHole(std::vector<int> const &src, int orb_pair, double &cum_sum, double &cpt);
    std::vector<int> selectDoubleHoles(std::vector<int> const &src,
    		std::vector<double> &cum_sum, std::vector<double> &cpt);

    // pick a holes for an opposite spin pair
    double oppSpinPairContribution(int i, int j, int a, int b);
    // contribution for a pair of same spin electrons
    double sameSpinPairContribution(int i, int j, int a, int b);

private:

    // The weighted selector knows the Hamiltonian
    Hamiltonian const &H;
    // and the starting point
    detType const &source;
    std::size_t norbs;

    // we do not copy/move WeightedSelectors
    WeightedSelector& operator=(WeightedSelector const&);
    WeightedSelector& operator=(WeightedSelector &&);
};

}

#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_WEIGHTEDSELECTOR_HPP_ */
