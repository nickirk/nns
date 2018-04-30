/*
 * SubspaceCF.hpp
 *
 *  Created on: Apr 10, 2018
 *      Author: guther
 */

#ifndef SRC_COSTFUNCTIONS_SUBSPACECF_HPP_
#define SRC_COSTFUNCTIONS_SUBSPACECF_HPP_

#include <Eigen/Dense>
#include "../Hamiltonian/SparseHMatrix.hpp"
#include "CostFunction.hpp"
#include "../Hamiltonian/Hamiltonian.hpp"
#include "../utilities/State.hpp"

namespace networkVMC{

// This cost function diagonalizes H in the subspace spanned by the
// sampled basis states and takes the distance of the resulting lowest
// eigenvector to the input state as cost function
class SubspaceCF: public CostFunction {
public:
	SubspaceCF(Hamiltonian const &H_):CostFunction(),H(H_),distance(0),subspaceEnergy(coeffType()){};
	virtual ~SubspaceCF();
// Derivative with respect to the input's coefficients
	virtual std::vector<Eigen::VectorXd > nabla(State const &input) const;
// Value of the cost function
	virtual double calc(State const &input) const {return distance;}
private:
	Hamiltonian const &H;
	double distance;
	coeffType subspaceEnergy;
// auxiliary function for getting the ground state in the space spanned by the determinants of input
	State diagonalizeSubspace(State const & input) const;
};

}

#endif /* SRC_COSTFUNCTIONS_SUBSPACECF_HPP_ */
