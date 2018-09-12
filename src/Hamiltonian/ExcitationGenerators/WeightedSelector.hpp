/*
 * WeightedSelector.hpp
 *
 *  Created on: Jun 19, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_EXCITATIONGENERATORS_WEIGHTEDSELECTOR_HPP_
#define SRC_HAMILTONIAN_EXCITATIONGENERATORS_WEIGHTEDSELECTOR_HPP_

#include "../../utilities/TypeDefine.hpp"


namespace networkVMC{

class TwoBodyHamiltonian;

/**
 * \class WeightedSelector
 * \brief Auxiliary excitation generation functionality
 *
 * Given a Hamiltonian and a determinant, we can select electrons based
 * on the Hamiltonian matrix elements.
*/

// it is essentially a wrapper to create cumulative sums
class WeightedSelector {
public:
	// this is created on the fly to select some parts/holes from a determinant
	// called source, weighted with a Hamiltonian H

	/**
	 * \param H_ Hamiltonian defining connectivity
	 * \param source_ basis vector to select from
	 */
	WeightedSelector(TwoBodyHamiltonian const &H_, detType const &source_):
		H(H_),source(source_),norbs(source_.size()){};
	virtual ~WeightedSelector();
    /**
     * \brief select a hole for a single excitation
     * \param[in] src electron to excite
     * \param[out] pgen probability of picking this hole
     * \param[out] hel matrix element of that excitation
     * \return index of a randomly chosen unoccupied orbital
     */
    int selectSingleHole(int src, double &pgen, double &hel);
    // select two holes based on a cumulative list
    int selectDoubleHole(std::vector<int> const &src, int orb_pair, double &cum_sum, double &cpt);
    std::vector<int> selectDoubleHoles(std::vector<int> const &src,
    		std::vector<double> &cum_sum, std::vector<double> &cpt);

    // pick a holes for an opposite spin pair
    double oppSpinPairContribution(int i, int j, int a, int b);
    // contribution for a pair of same spin electrons
    double sameSpinPairContribution(int i, int j, int a, int b);

    /// assignment is deleted
    WeightedSelector& operator=(WeightedSelector const&) = delete;
    /// assignment is deleted
    WeightedSelector& operator=(WeightedSelector &&) = delete;
private:

    /// The underlying Hamiltonian
    TwoBodyHamiltonian const &H;
    /// The starting point of selection
    detType const &source;
    /// Number of single-particle basis functions (orbitals)
    std::size_t norbs;
};

}

#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_WEIGHTEDSELECTOR_HPP_ */
