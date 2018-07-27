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
#include <vector>
#include <set>
#include <memory>

namespace networkVMC {

// generic Hamiltonian on a lattice

class LatticeHamiltonian: public Hamiltonian {
public:
	// pass the dimensions of the lattice as variadic template arguments
	// currently only square lattices are supported, but this might change
	template<typename ...Args>
	// construct the lattice from args
	LatticeHamiltonian(Args ...args):grid(std::unique_ptr<Lattice>(new Lattice(std::forward<Args>(args)...) ) ){};
	virtual ~LatticeHamiltonian(){};

	// number of sites in grid
	int size() const {return grid->size();}
	// adjacent sites (for testing purposes)
	std::set<int> adjacents(int i){return grid->adjacents(i);}

	// generate coupled states on a lattice
	virtual std::vector<detType> getCoupledStates(detType const &source) const;
protected:
	// the lattice on which the Hamiltonian is defined
	// use a ptr to make introduction of polymorphic lattice easy
	std::unique_ptr<Lattice> grid;
};

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_LATTICEHAMILTONIAN_HPP_ */
