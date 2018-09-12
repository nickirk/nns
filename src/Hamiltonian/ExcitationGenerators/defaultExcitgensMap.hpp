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

/// Pointer to ExcitationGenerator
using excitgenPtr = std::unique_ptr<ExcitationGenerator>;

/**
 * \fn std::unique_ptr<ExcitationGenerator> getDefaultExcitgen(Hamiltonian const &H, detType const &HF);
 * \brief Return the default ExcitationGenerator for a given Hamiltonian
 * \param[in] H Hamiltonian for which to get the ExcitationGenerator
 * \param[in] HF reference basis vector for the ExcitationGenerator (required for initialization)
 * \return unique pointer to a new ExcitationGenerator of the type corresponding to H
 */
excitgenPtr getDefaultExcitgen(Hamiltonian const &H,
		detType const &HF);

}

#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_DEFAULTEXCITGENSMAP_HPP_ */
