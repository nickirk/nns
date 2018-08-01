/*
 * LocalLatticeHamiltonian.hpp
 *
 *  Created on: Jul 27, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_LATTICE_HPP_
#define SRC_HAMILTONIAN_LATTICE_HPP_

#include <set>
#include <vector>
#include <numeric>

namespace networkVMC {

// lattice class for LatticeHamiltonians
// the lattice has no idea about the site-local hilbert space

class Lattice {
public:
	// variadic template allows to construct with arbitrary spatial dimensions supplied
	// initializer_list would do the same, but forces to use {}-initialization
	// maybe switch if not needed
	template<typename ...Args>
	Lattice(bool pbc_, Args... args):latticeDim({args...}),pbc(pbc_){
		// if no template arguments are supplied, we assume a 1-d model (as that is the only one where we
		// can build the hamiltonian without knowing the dimensions)
		if(latticeDim.size()==0) latticeDim.push_back(12);
		// number of sites in total
		numSites = std::accumulate(latticeDim.begin(),latticeDim.end(),1,std::multiplies<int>());
		for(int i = 0; i < numSites; ++i){
			adjacencyList.push_back( adjacentSites(i) );
		}
	};

	// access utility for usage in LatticeHamiltonian:
	// get the number of spatial sites in the lattice
	int size() const {return numSites;}
	// get the spatial dimension of the lattice
	int dimension() const {return latticeDim.size();}
	// get the sites adjacent to site i
	std::set<int>const & adjacents(int i) const {return adjacencyList[i];}
	// if two sites with indices i, j are adjacent
	bool isAdjacent(int i, int j) const;

	// we might want to generalize to other lattices than square at some point
	virtual ~Lattice(){};
private:
	// spatial dimensions of the lattice
	std::vector<int> latticeDim;
	int numSites;
	// all sites adjacent to i
	std::set<int> adjacentSites(int i) const;
	// adjacent sites for each sites, for faster runtime evaluation of matrix elements
	std::vector< std::set<int> > adjacencyList;
	// if the lattice has periodic boundary
	bool pbc;
};


} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_LATTICE_HPP_ */
