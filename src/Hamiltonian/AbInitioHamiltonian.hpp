/*
 * AbInitioHamiltonian.h
 *
 *  Created on: Feb 13, 2018
 *      Author: Lauretta Schwarz
 */

#ifndef SRC_HAMILTONIAN_ABINITIOHAMILTONIAN_HPP_
#define SRC_HAMILTONIAN_ABINITIOHAMILTONIAN_HPP_

#include <string>
#include "Hamiltonian.hpp"

namespace networkVMC{

class AbInitioHamiltonian: public Hamiltonian {
public:
	AbInitioHamiltonian(int dimension):Hamiltonian(dimension){};
	virtual ~AbInitioHamiltonian();
	int getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex) const;
    // count the number of connected states
    int countNumberCoupledStates(detType const &source, int exflag, int &nsingleexcit, int &ndoubleexcit);
    // calculate the generation probability
    double calcGenProp(detType const &source, detType const &target);
    // return the next electron pair (for the deterministic excitation generator)
	void getElecPair(std::vector<int> const &source_orbs, int &el1, int &el2, int &spinpair, int &ind) const;
    // return the next single excitation (deterministic excitation generator)
    detType getSingleExcitation(detType const &source, std::vector<int> const &source_orbs, std::vector<int> &holes, std::vector<int> &particles, int &exflag, bool &ballexcitfound) const;
    // return the next double excitation (deterministic excitation generator)
    detType getDoubleExcitation(detType const &source, std::vector<int> const &source_orbs, std::vector<int> &holes, std::vector<int> &particles, int &exflag, bool &ballexcitfound) const;
    // pick a random electron pair (random excitation generator)
    std::vector<int> pickElecPair(std::vector<int> const &source_orbs, std::vector<int> &spin, int &elecpairs) const;
    // could the number of forbidden spin orbitals (random excitation generator)
    int countForbiddenOrbs(std::vector<int> const &spins, std::vector<int> nunoccs) const;
    // pick orbital a of an ij->ab excitation (random excitation generator)
    int pickOrbA(detType const &source, std::vector<int> const &spins, std::vector<int> noccs, std::vector<int> unoccs, int norbs, int nel, int nforbiddenorbs, int &nexcita, int &aspin, bool &baorbfail) const;
    // pick orbital b of an ij->ab excitation (random excitation generator)
    int pickOrbB(detType const &source, std::vector<int> const &spins, std::vector<int> noccs, std::vector<int> unoccs, int norbs, int nel, int nforbiddenorbs, int aorb, int aspin, int &nexcitb, int &bspin, int &nexcitotherway) const;
    // randomly generate a single excitation excitation (random excitation generator)
    detType genSingleExcitation(detType const &source, std::vector<int> const &source_orbs, std::vector<int> &holes, std::vector<int> &particles, double &pGet, double pdoubnew, std::vector<int> noccs, std::vector<int> nunoccs) const;
    // randomly generate a double excitation excitation (random excitation generator)
    detType genDoubleExcitation(detType const &source, std::vector<int> const &source_orbs, std::vector<int> &holes, std::vector<int> &particles, double &pGet, double pdoubnew, std::vector<int> noccs, std::vector<int> nunoccs) const;
    // random excitation generator: randomly generates one connected determinant
    detType getRandomCoupledState(detType const &source, double &p) const;
    // deterministic excitation generator: generates all connected determinants
    std::vector<detType> getCoupledStates(detType const &source) const;
};

AbInitioHamiltonian readAbInitioHamiltonian(int dim, std::string file_name);

}

#endif /* SRC_HAMILTONIAN_ABINITIOHAMILTONIAN_HPP_ */
