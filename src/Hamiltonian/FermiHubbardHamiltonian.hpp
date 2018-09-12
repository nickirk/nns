/*
 * FermiHubbardHamiltonian.hpp
 *
 *  Created on: Jun 20, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_FERMIHUBBARDHAMILTONIAN_HPP_
#define SRC_HAMILTONIAN_FERMIHUBBARDHAMILTONIAN_HPP_

#include "LatticeHamiltonian.hpp"

namespace networkVMC {

/**
 * \class FermiHubbardHamiltonian
 * \brief Implements the Fermi-Hubbard Hamiltonian
 *
 * Implementation of the fermionic Hubbard model, based on the LatticeHamiltonian. Contains a local on-site interaction
 * and a nearest-neighbor hopping.
 */
class FermiHubbardHamiltonian: public LatticeHamiltonian {
  public:

	// generate an n-d lattice
	/**
	 * \tparam ...Args types of the lattice dimensions
	 * \param[in] U_ on-site Hubbard interaction
	 * \param[in] t_ hopping parameter
	 * \param[in] ...args lattice dimensions
	 */
	template<typename ...Args>
	FermiHubbardHamiltonian(double U_, double t_, Args ...args):
	LatticeHamiltonian(true,args...),U(U_),t(t_){};

	// get the hubbard matrix element
	double operator()(detType const &a, detType const &b) const;
	virtual ~FermiHubbardHamiltonian(){};

	/// \return Hubbard
	virtual HType type() const {return Hubbard;}
private:
	/// Parameters of the matrix elements
	double U, t;
	void addCoupledStates(std::vector<detType> &list, detType const &source, int siteA, int siteB) const;
};

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_FERMIHUBBARDHAMILTONIAN_HPP_ */
