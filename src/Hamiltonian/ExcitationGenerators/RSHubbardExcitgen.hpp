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

class RSHubbardExcitgen: public clonableExcitgen<RSHubbardExcitgen> {
public:
	RSHubbardExcitgen();
	virtual ~RSHubbardExcitgen();
	virtual detType generateExcitation(
			detType const &source,  double &pGet);
	virtual double getExcitationProb(
			detType const &source, detType const &target);
};

void getRSHubSpawnLists(detType const &source, std::vector<int> &spawnLeft,
		std::vector<int> &spawnRight);

} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_RSHUBBARDEXCITGEN_HPP_ */
