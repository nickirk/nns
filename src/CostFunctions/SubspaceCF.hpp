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

namespace networkVMC{

// This cost function diagonalizes H in the subspace spanned by the
// sampled basis states and takes the distance of the resulting lowest
// eigenvector to the input state as cost function
class SubspaceCF: public CostFunction {
public:
	SubspaceCF(Hamiltonian const &H_);
	virtual ~SubspaceCF();
// Derivative with respect to the input's coefficients
	virtual std::vector<Eigen::VectorXd > nabla(std::vector<State> const &input) const;
// Value of the cost function
	virtual double calc(std::vector<State> const &input) const;
private:
	Hamiltonian const &H;
	SparseHMatrix HS;
// auxiliary function for getting the ground state in the space spanned by the determinants of input
	std::vector<State> diagonalizeSubspace(std::vector<State> const & input);
};

}

#endif /* SRC_COSTFUNCTIONS_SUBSPACECF_HPP_ */
