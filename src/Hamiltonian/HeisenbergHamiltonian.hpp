/*
 * HeisenbergHamiltonian.hpp
 *
 *  Created on: Jul 26, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_HEISENBERGHAMILTONIAN_HPP_
#define SRC_HAMILTONIAN_HEISENBERGHAMILTONIAN_HPP_

#include "TwoBodyHamiltonian.hpp"
#include <set>
#include <numeric>

namespace networkVMC {

class HeisenbergHamiltonian: public Hamiltonian {
public:
	// variadic template allows to construct with arbitrary spatial dimensions supplied
	// initializer_list would do the same, but forces to use {}-initialization
	// maybe switch if not needed
	template<typename ...Args>
	HeisenbergHamiltonian(double J_, Args... args):
		Hamiltonian(),latticeDim({args...}),J(J_){
		// if no template arguments are supplied, we assume a 1-d model (as that is the only one where we
		// can build the hamiltonian without knowing the dimensions)
		if(latticeDim.size()==0) latticeDim.push_back(12);
		// number of sites in total
		numSites = std::accumulate(latticeDim.begin(),latticeDim.end(),1,std::multiplies<int>());
		for(int i = 0; i < numSites; ++i){
			adjacencyList.push_back( adjacentSites(i) );
		}
	};

	// get the matrix element between a and b
	double operator()(detType const &a, detType const &b) const;
	// generate all states coupled to source
	std::vector<detType> getCoupledStates(detType const &source) const;

	virtual ~HeisenbergHamiltonian(){};
	HType type() const {return Heisenberg;}

	int size() const {return numSites;}
	std::set<int> adj(int i) const {return adjacencyList[i];}
private:
	// spatial dimensions of the lattice
	std::vector<int> latticeDim;
	int numSites;
	// coupling parameter
	double J;
	// we cannot set matrix elements for the Heisenberg Hamiltonian
	void setMatrixElement(int a, int b, double newEntry);
	void setMatrixElement(int a, int b, int c, int d, double newEntry);
	// if two sites with indices i, j are adjacent
	bool isAdjacent(int i, int j) const;
	// all sites adjacent to i
	std::set<int> adjacentSites(int i) const;
	// adjacent sites for each sites, for faster runtime evaluation of matrix elements
	std::vector< std::set<int> > adjacencyList;
};

//---------------------------------------------------------------------------------------------------//

// convert an occupancy flag to the ms value (0 is down, 1 is up)
inline int msVal(bool spin){
	if(spin) return 1;
	return -1;
}

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_HEISENBERGHAMILTONIAN_HPP_ */
