/*
 * Hamiltonian.hpp
 *
 *  Created on: Jul 26, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_HAMILTONIAN_HPP_
#define SRC_HAMILTONIAN_HAMILTONIAN_HPP_

#include "../utilities/TypeDefine.hpp"

namespace networkVMC {

// generic base class for Hamiltonians

/**
 * \class Hamiltonian
 * \brief Abstract base class for Hamiltonians
 *
 * This class provides the interface for implementing Hamiltonians. A Hamiltonian is defined by its operator()
 * which returns the matrix element between two basis vectors.
 * An additionally required function is the getCoupledStates member function that returns all basis vectors that
 * have nonzero matrix element with an input basis vector
 * Keeping it general, Hamiltonians do not have information on the dimension of the space from which the
 * basis vectors come, as their implementation might not actually depend on that.
 */
class Hamiltonian {
  public:
	Hamiltonian();
	/**
	 * \brief Gets the matrix element between two basis states
	 *
	 * \param a first basis state
	 * \param b second basis state
	 * \return matrix element of the Hamiltonian between a and b
	 * Hamiltonians need to be able to produce a matrix element given two determinants
	 */
	virtual double operator()(detType const &a, detType const &b) const =0;
	// remarkably, what they do not have, is information on the size of the input determinants
	// in general, Hamiltonians can treat arbitrary a,b, but specifications might constrain this
	virtual ~Hamiltonian();

	/**
	 * \brief Generate all states coupling to source
	 *
	 * \param[in] source basis vector to which to couple
	 * \return list of all basis vectors coupling to source
	 * Note that - in contrast to excitation generation - this is really
	 * a property of the Hamiltonian
	 */
	virtual std::vector<detType> getCoupledStates(detType const &source) const = 0;

	/**
	 * \brief type of the Hamiltonian
	 *
	 * \return enum indicating of which type the Hamiltonian is
	 *
	 * Possible types are Constant, Hubbard, AbInitio and Heisenberg. Default is Constant
	 */
	virtual HType type() const{return Constant;};
};

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_HAMILTONIAN_HPP_ */
