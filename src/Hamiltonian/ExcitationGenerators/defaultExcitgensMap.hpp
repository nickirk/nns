/*
 * defaultExcitgensMap.hpp
 *
 *  Created on: Jun 21, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_EXCITATIONGENERATORS_DEFAULTEXCITGENSMAP_HPP_
#define SRC_HAMILTONIAN_EXCITATIONGENERATORS_DEFAULTEXCITGENSMAP_HPP_

#include "../../utilities/TypeDefine.hpp"
#include <memory>

namespace networkVMC{

class ExcitationGenerator;
class Hamiltonian;

using excitgenPtr = std::unique_ptr<ExcitationGenerator>;

excitgenPtr getDefaultExcitgen(Hamiltonian const &H,
		detType const &HF);

}

#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_DEFAULTEXCITGENSMAP_HPP_ */
