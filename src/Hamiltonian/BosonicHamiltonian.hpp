/*
 * BosonicHamiltonian.hpp
 *
 *  Created on: Nov 15, 2017
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_BOSONICHAMILTONIAN_HPP_
#define SRC_HAMILTONIAN_BOSONICHAMILTONIAN_HPP_

#include "TwoBodyHamiltonian.hpp"

namespace networkVMC{

/**
 * \class BosonicHamiltonian
 * \brief TwoBodyHamiltonian with bosoinc commutation relations
 */
class BosonicHamiltonian: public TwoBodyHamiltonian {
  public:
	/// \params[in] dimension dimension of the one-particle Hilbert space
	BosonicHamiltonian(int dimension):TwoBodyHamiltonian(dimension){};
	virtual ~BosonicHamiltonian();
/// \return 1
	int getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex) const{return 0;};
    //detType getRandomCoupledState(detType const &source, double &p){};
    std::vector<detType> getCoupledStates(detType const &source) const;
};

}

#endif /* SRC_HAMILTONIAN_BOSONICHAMILTONIAN_HPP_ */
