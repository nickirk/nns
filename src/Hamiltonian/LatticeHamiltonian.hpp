/*
 * LatticeHamiltonian.hpp
 *
 *  Created on: Jul 27, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_LATTICEHAMILTONIAN_HPP_
#define SRC_HAMILTONIAN_LATTICEHAMILTONIAN_HPP_

#include "Hamiltonian.hpp"
#include "Lattice.hpp"
#include "../utilities/TypeDefine.hpp"
#include "../utilities/DeepCpyUniquePtr.hpp"
#include <vector>
#include <set>
#include <memory>

namespace networkVMC {

/**
 * \class LatticeHamiltonian
 * \brief Generic Hamiltonian on a lattice
 *
 * A Hamiltonian whose couplings are not fully stored as in the TwoBodyHamiltonian, but
 * have a structure given by a Lattice, which defines which single-particle states interact with which other ones.
 */

class LatticeHamiltonian: public Hamiltonian {
  public:
	// pass the dimensions of the lattice as variadic template arguments
	// currently only square lattices are supported, but this might change
	/**
	 * \tparam ...Args types of the lattice dimensions
	 * \param[in] ...args dimension of the lattice
	 */
	template<typename ...Args>
	// construct the lattice from args
	LatticeHamiltonian(Args ...args):grid(DeepCpyUniquePtr<Lattice>(
			new Lattice(std::forward<Args>(args)...) ) ){};
	virtual ~LatticeHamiltonian(){};

	/// \return number of sites in grid
	int size() const {return grid->size();}
	/// \return spatial dimension of the grid
	int dimension() const {return grid->dimension();}
	/**
	 * \brief Get the sites adjacent to site i (for testing purposes)
	 * \param[in] i Site to consider
	 * \return List of adjacent sites (by lattice index)
	 */
	std::set<int>const & adjacents(int i) const{return grid->adjacents(i);}

	// generate coupled states on a lattice
	virtual std::vector<detType> getCoupledStates(detType const &source) const;
protected:

	// use a ptr to make introduction of polymorphic lattice easy
	/// The lattice on which the Hamiltonian is defined
	DeepCpyUniquePtr<Lattice> grid;
	// we re-use this in derived classes to get the coupled states, as the derived versions
	// boil down to repeated calls to this one

	/**
	 * \brief add certain basis vectors to a list
	 * \param list list to append to
	 * \param[in] source basis vector to refer to
	 * \param[in] siteA, siteB sites to consider
	 * Adds all states coupling to source via siteA and siteB to list
	 */
	virtual void addCoupledStates(std::vector<detType> &list, detType const &source, int siteA, int siteB) const;
};

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_LATTICEHAMILTONIAN_HPP_ */
