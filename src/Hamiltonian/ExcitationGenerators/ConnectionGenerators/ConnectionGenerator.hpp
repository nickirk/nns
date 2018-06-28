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

// this takes a determinant, a Hamiltonian and a number and generates a list of that many connected states
// pGen is the array of probabilities
std::vector<detType> sampleConnections(Hamiltonian const &H, detType const &base,
		int numConnections, std::vector<double> &pGen);


} /* namespace networkVMC */

#endif /* SRC_HAMILTONIAN_EXCITATIONGENERATORS_CONNECTIONGENERATORS_CONNECTIONGENERATOR_HPP_ */
