/*
 * AbInitioHamiltonian.h
 *
 *  Created on: Feb 13, 2018
 *      Author: Lauretta Schwarz
 */

#ifndef SRC_HAMILTONIAN_ABINITIOHAMILTONIAN_HPP_
#define SRC_HAMILTONIAN_ABINITIOHAMILTONIAN_HPP_

#include <string>
#include "FermionicHamiltonian.hpp"

namespace networkVMC{

class AbInitioHamiltonian: public FermionicHamiltonian {
public:
    AbInitioHamiltonian(int dimension):FermionicHamiltonian(dimension){};
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
            int &el1, int &el2, int &spinpair, int &ind) const;

    // deterministic excitation generator: generates all connected determinants
    std::vector<detType> getCoupledStates(detType const &source) const;
    // calculate the generation probability biased according to the connecting
    // Hamiltonian matrix element

    // typecheck for setting defaults
    virtual HType type() const {return AbInitio;}
};

AbInitioHamiltonian readAbInitioHamiltonian(int dim, std::string file_name);

}

#endif /* SRC_HAMILTONIAN_ABINITIOHAMILTONIAN_HPP_ */
