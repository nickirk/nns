/*
 * AbInitioHamiltonian.h
 *
 *  Created on: Feb 13, 2018
 *      Author: Lauretta Schwarz
 */

#ifndef SRC_HAMILTONIAN_ABINITIOHAMILTONIAN_HPP_
#define SRC_HAMILTONIAN_ABINITIOHAMILTONIAN_HPP_

#include <string>
#include "TwoBodyHamiltonian.hpp"
#include "../HilbertSpace/Determinant.hpp"

namespace networkVMC{

// a fermionic generic two-body hamiltonian (aka Ab-initio Hamiltonian)
class AbInitioHamiltonian: public TwoBodyHamiltonian {
  public:
	// dimension being the number of orbitals
    AbInitioHamiltonian(int dimension):TwoBodyHamiltonian(dimension){};
	virtual ~AbInitioHamiltonian();
    // count the number of connected states
    int countNumberCoupledStates(detType const &source, int exflag, int &nsingleexcit, int &ndoubleexcit);

    // deterministic excitgen
    // TODO: move to its own class
    detType getSingleExcitation(detType const &source,
    		std::vector<int> const &source_orbs, std::vector<int> &holes,
			std::vector<int> &particles, int &exflag, bool &ballexcitfound) const;

    detType getDoubleExcitation(detType const &source, std::vector<int>
    const &source_orbs, std::vector<int> &holes, std::vector<int> &particles,
	int &exflag, bool &ballexcitfound) const;

    void getElecPair(std::vector<int> const &source_orbs,
            int &el1, int &el2, int &spinpair, int ind) const;

    // deterministic excitation generator: generates all connected determinants
    std::vector<detType> getCoupledStates(detType const &source) const;
    // calculate the generation probability biased according to the connecting
    // Hamiltonian matrix element

    // typecheck for setting defaults
    virtual HType type() const {return AbInitio;}

    // The fermi sign is the main difference to the plain TwoBodyHamiltonian
    int getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex) const{
    	return JWStringLength(alpha,annihilatorIndex,creatorIndex);
    }
};

AbInitioHamiltonian readAbInitioHamiltonian(std::string file_name, bool molpro_fcidump=false);

}

#endif /* SRC_HAMILTONIAN_ABINITIOHAMILTONIAN_HPP_ */
