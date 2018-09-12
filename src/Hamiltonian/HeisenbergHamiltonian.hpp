/*
 * HeisenbergHamiltonian.hpp
 *
 *  Created on: Jul 26, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_HEISENBERGHAMILTONIAN_HPP_
#define SRC_HAMILTONIAN_HEISENBERGHAMILTONIAN_HPP_

#include "LatticeHamiltonian.hpp"

namespace networkVMC {

/**
 * \class HeisenbergHamiltonian
 * \brief Implements a Heisenberg spin-1/2 Hamiltonian
 *
 * Lattice-based spin-1/2 hamiltonian with isotropic nearest-neighbor interaction,
 * commonly known as spin-1/2 Heisenberg Hamiltonian
 */
class HeisenbergHamiltonian: public LatticeHamiltonian {
  public:
	/**
	 * \tparam ...Args Types of the dimensions of the lattice
	 * \param[in] J_ Value of the coupling constant (determines sign of the interaction)
	 * \param[in] pbc Flag to indicate periodic boundaries
	 * \param[in] ...args Dimensions of the lattice
	 */
	template<typename ...Args>
	HeisenbergHamiltonian(double J_, bool pbc, Args ...args):
		LatticeHamiltonian(pbc, args...),J(J_){};
	/**
	 * \overload
	 * Sets the pbc flag to true
	 */
	template<typename ...Args>
	HeisenbergHamiltonian(double J_, Args ...args):
		LatticeHamiltonian(true, args...),J(J_){};


	double operator()(detType const &a, detType const &b) const;

	/// \return Heisenberg
	HType type() const {return Heisenberg;}
	virtual ~HeisenbergHamiltonian(){};
private:
	double J;
};

//---------------------------------------------------------------------------------------------------//

// convert an occupancy flag to the ms value (0 is down, 1 is up)
/**
 * \fn inline int msVal(bool spin)
 * \param[in] spin 'occupation' value from a basis vector
 * \return Sign of ms of that state when converted to spin-1/2
 */
inline int msVal(bool spin){
	if(spin) return 1;
	return -1;
}

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_HEISENBERGHAMILTONIAN_HPP_ */
