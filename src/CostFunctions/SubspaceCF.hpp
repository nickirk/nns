/*
 * SubspaceCF.hpp
 *
 *  Created on: Apr 10, 2018
 *      Author: guther
 */

#ifndef SRC_COSTFUNCTIONS_SUBSPACECF_HPP_
#define SRC_COSTFUNCTIONS_SUBSPACECF_HPP_

#include <Eigen/Dense>
#include "CostFunction.hpp"
#include "../utilities/State.hpp"

namespace networkVMC{

class Hamiltonian;

// This cost function diagonalizes H in the subspace spanned by the
// sampled basis states and takes the distance of the resulting lowest
// eigenvector to the input state as cost function
class SubspaceCF: public CostFunction {
public:
	// SubspaceCF creation/destruction
	SubspaceCF(Hamiltonian const &H_):CostFunction(),H(H_),distance(0),subspaceEnergy(coeffType()){};
	virtual ~SubspaceCF();
// Derivative with respect to the input's coefficients
	nablaType nabla(State const &input) const;
// Value of the cost function
	double calc(State const &input) const {return distance;}
private:
  // The underlying Hamiltonian
	Hamiltonian const &H;
  // cache variable for the distance between input state and subspace eigenstate
	mutable double distance;
  // cache variable for subspace energy
	mutable coeffType subspaceEnergy;
  // auxiliary function for getting the ground state in the space spanned by the determinants of input
	State diagonalizeSubspace(State const & input) const;

  // Has a reference member, so assignment is not a thing
  SubspaceCF& operator=(SubspaceCF const &source);
  SubspaceCF& operator=(SubspaceCF &&source);
};

}

#endif /* SRC_COSTFUNCTIONS_SUBSPACECF_HPP_ */
