/*
 * RSHubbardExcitgen.hpp
 *
 *  Created on: Jun 19, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_EXCITATIONGENERATORS_RSHUBBARDEXCITGEN_HPP_
#define SRC_HAMILTONIAN_EXCITATIONGENERATORS_RSHUBBARDEXCITGEN_HPP_

#include "ExcitationGenerator.hpp"
#include "../../HilbertSpace/Determinant.hpp"

namespace networkVMC {

/**
 * \class RSHubbardExcitgen
 * \brief ExcitationGenerator specifically for 1-d real-space Hubbard systems
 */
class RSHubbardExcitgen: public ClonableExcitgen<RSHubbardExcitgen> {
public:
	RSHubbardExcitgen();
	virtual ~RSHubbardExcitgen();
	detType generateExcitation(
			detType const &source,  double &pGet);
	double getExcitationProb(
			detType const &source, detType const &target);
};

/**
 * \fn void getRSHubSpawnLists(detType const &source, std::vector<int> &spawnLeft, std::vector<int> &spawnRight);
 * \brief The RSHubbardExcitgen version of the Lattice member function adjacents
 * \param[in] source basis vector to spawn from
 * \param[out] spawnLeft list of sites that can be spawned upon in a left-move
 * \param[out] spawnRight list of sites that can be spawned upon in a right-move
 *
 * For a given basis vector, creates a list of all sites that can be reached by moving an electro to
 * the left/right.
 */
void getRSHubSpawnLists(detType const &source, std::vector<int> &spawnLeft,
		std::vector<int> &spawnRight);

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_RSHUBBARDEXCITGEN_HPP_ */
