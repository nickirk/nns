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

/**
 * \class Lattice
 * \brief Lattice class for LatticeHamiltonians
 */
// the lattice has no idea about the site-local hilbert space

class Lattice {
  public:
	// variadic template allows to construct with arbitrary spatial dimensions supplied
	// initializer_list would do the same, but forces to use {}-initialization
	// maybe switch if not needed

	/**
	 * \tparam ...Args Types of the lattice dimension
	 * \param[in] pbc_ Flag to indicate periodic boundary conditions
	 * \param[in] ...args Dimensions of the lattice
	 */
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
	/// \return The number of spatial sites in the lattice
	int size() const {return numSites;}
	/// \return The spatial dimension of the lattice
	int dimension() const {return latticeDim.size();}

	/**
	 * \brief Get the sites adjacent to site i
	 * \param[in] i Site to consider
	 * \return List of adjacent sites (by lattice index)
	 * This function gets the adjacent sites from a cache
	 */
	std::set<int>const & adjacents(int i) const {return adjacencyList[i];}
	/**
	 * \brief Checks if two sites with indices i, j are adjacent
	 * \param[in] i,j lattice indices of the sites
	 * \return Flag indicating if i,j are adjacent
	 */
	bool isAdjacent(int i, int j) const;

	// we might want to generalize to other lattices than square at some point
	virtual ~Lattice(){};
	/// Dynamic polymorphic copy
	virtual Lattice* clone() const {return new Lattice(*this);}
private:
	/// Spatial dimensions of the lattice
	std::vector<int> latticeDim;
	/// Number of lattice sites
	int numSites;
	/**
	 * \brief Generate a list of adjacent sites from scratch
	 * \param[in] i Site to consider
	 * \return List of adjacent sites (by lattice index)
	 * This function generates the list of adjacent sites from scratch,
	 * it is used to create the cache used in adjacents()
	 */
	std::set<int> adjacentSites(int i) const;
	/// Adjacent sites for each sites, for faster runtime evaluation of matrix elements
	std::vector< std::set<int> > adjacencyList;
	/// Flag indicating if the lattice has periodic boundary
	bool pbc;
};


} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_LATTICE_HPP_ */
