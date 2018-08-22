/*
 * inputStateGenerator.hpp
 *
 *  Created on: Jun 28, 2018
 *      Author: guther
 */

#ifndef SRC_UTILITIES_INPUTSTATEGENERATOR_HPP_
#define SRC_UTILITIES_INPUTSTATEGENERATOR_HPP_

#include "State.hpp"
#include "../Network/Parametrization.hpp"
#include "../Samplers/Sampler.hpp"

namespace networkVMC {

// forward declarations
//class Sampler;
class Hamiltonian;

//class Parametrization;

// this class takes a sampler and a hamiltonian and creates the sampled states
// (with or without connections)

template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class InputStateGenerator {
  public:
	InputStateGenerator(Sampler<coeffType> &msampler_, Hamiltonian const &H_, Parametrization<F, coeffType> const &para);
	virtual ~InputStateGenerator();

	// create an input state, with connections indicating if connections are to be sampled
	State<coeffType> generate(int numCons) const;
private:
	// the sampler used for getting random states
	Sampler<coeffType> &msampler;
	// the Hamiltonian governing the connections
	Hamiltonian const &H;
	// the parametrization containing the coefficients of the determinants
	Parametrization<F, coeffType> const &para;

};

} /* namespace networkVMC */

#endif /* SRC_UTILITIES_INPUTSTATEGENERATOR_HPP_ */
