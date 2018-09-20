/*
 * inputStateGenerator.hpp
 *
 *  Created on: Jun 28, 2018
 *      Author: guther
 */

#ifndef SRC_UTILITIES_INPUTSTATEGENERATOR_HPP_
#define SRC_UTILITIES_INPUTSTATEGENERATOR_HPP_

#include "State.hpp"
#include "../Samplers/Sampler.hpp"

namespace networkVMC {


class Hamiltonian;
class Parametrization;

/**
 * \class InputStateGenerator
 * \brief Creates randomly sampled States
 *
 * This class takes a Sampler and a Hamiltonian and creates the sampled states
 * (with or without connections)
 */
class InputStateGenerator {
  public:
	/**
	 * \param msampler_ Sampler to use for basis vector selection
	 * \param[in] H_ Hamiltonian defining connectivity
	 * \param[in] para Parametrization of the vector coefficients
	 */
	InputStateGenerator(Sampler &msampler_, Hamiltonian const &H_, Parametrization const &para);
	virtual ~InputStateGenerator();

	/*
	 * \brief create an input state
	 * \param[in] numCons number of connected basis vectors to include per sampled basis vector
	 * \return State containing the sampled basis vectors and their coefficients
	 */
	State generate(int numCons) const;
private:
	/// the Sampler used for getting random states
	Sampler &msampler;
	/// the Hamiltonian governing the connections
	Hamiltonian const &H;
	/// the Parametrization containing the coefficients of the determinants
	Parametrization const &para;

};

} /* namespace networkVMC */

#endif /* SRC_UTILITIES_INPUTSTATEGENERATOR_HPP_ */
