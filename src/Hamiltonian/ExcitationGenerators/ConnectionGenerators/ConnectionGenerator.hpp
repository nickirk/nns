/*
 * ConnectionGenerator.hpp
 *
 *  Created on: Jun 28, 2018
 *      Author: guther
 */

#ifndef SRC_HAMILTONIAN_EXCITATIONGENERATORS_CONNECTIONGENERATORS_CONNECTIONGENERATOR_HPP_
#define SRC_HAMILTONIAN_EXCITATIONGENERATORS_CONNECTIONGENERATORS_CONNECTIONGENERATOR_HPP_

#include "../../../utilities/TypeDefine.hpp"
#include <vector>

namespace networkVMC {

class Hamiltonian;

//
//
/**
 * \brief Randomly generate a list of basis vectors connected to a given basis vector
 * \fn std::vector<detType> sampleConnections(Hamiltonian const &H, detType const &base, int numConnections, std::vector<double> &pGen);
 * \param[in] H hamiltonian defining connectivity
 * \param[in] base basis vector to sample from
 * \param[in] numConnections number of connected states to sample
 * \param[out] pGen array of generation probabilities
 * \return list of randomly generated basis vectors
 * Takes a determinant, a Hamiltonian and a number and generates a list of that many connected states
 */
std::vector<detType> sampleConnections(Hamiltonian const &H, detType const &base,
		int numConnections, std::vector<double> &pGen);


} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_CONNECTIONGENERATORS_CONNECTIONGENERATOR_HPP_ */
