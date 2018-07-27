/*
 * FermionBasis.hpp
 *
 *  Created on: Jul 26, 2018
 *      Author: guther
 */

#ifndef SRC_HILBERTSPACE_FERMIONBASIS_HPP_
#define SRC_HILBERTSPACE_FERMIONBASIS_HPP_

#include "../utilities/SpinConfig.hpp"
#include <vector>
#include "../utilities/TypeDefine.hpp"

namespace networkVMC {

std::vector<detType> generateFermionBasis(SpinConfig const &spinConfig);
void createFermionBasisDet(std::vector<detType> &basis, std::vector<int> &combination, int offset,
		int numEle, std::vector<int> const &listOfOrbNum);

} /* namespace networkVMC */

#endif /* SRC_HILBERTSPACE_FERMIONBASIS_HPP_ */
