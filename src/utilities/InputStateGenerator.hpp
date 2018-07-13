/*
 * inputStateGenerator.hpp
 *
 *  Created on: Jun 28, 2018
 *      Author: guther
 */

#ifndef SRC_UTILITIES_INPUTSTATEGENERATOR_HPP_
#define SRC_UTILITIES_INPUTSTATEGENERATOR_HPP_

#include "State.hpp"

namespace networkVMC {

// forward declarations
class Sampler;
class Hamiltonian;

template<typename T>
class Parametrization;

// this class takes a sampler and a hamiltonian and creates the sampled states
// (with or without connections)

template<typename T=VecType>
class InputStateGenerator {
public:
	InputStateGenerator(Sampler &msampler_, Hamiltonian const &H_, Parametrization<T> const &para);
	virtual ~InputStateGenerator();

	// create an input state, with connections indicating if connections are to be sampled
	State generate(int numCons) const;
private:
	// the sampler used for getting random states
	Sampler &msampler;
	// the Hamiltonian governing the connections
	Hamiltonian const &H;
	// the parametrization containing the coefficients of the determinants
	Parametrization<T> const &para;

};

} /* namespace networkVMC */

#endif /* SRC_UTILITIES_INPUTSTATEGENERATOR_HPP_ */
