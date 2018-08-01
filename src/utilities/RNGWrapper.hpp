/*
 * RNGWrapper.hpp
 *
 *  Created on: Jul 30, 2018
 *      Author: guther
 */

#ifndef SRC_UTILITIES_RNGWRAPPER_HPP_
#define SRC_UTILITIES_RNGWRAPPER_HPP_

#include <random>

namespace networkVMC {

// wrapper class for random number generation
class RNGWrapper {
public:
	RNGWrapper():rd(),rng(std::mt19937(rd())),uni(0,1){};
	virtual ~RNGWrapper(){};
	double operator()() const{return uni(rng);};
private:
	std::random_device rd;     // only used once to initialise (seed) engine
	mutable std::mt19937 rng;    // random-number engine used (Mersenne-Twister in this case)
	mutable std::uniform_real_distribution<double> uni;		// Uniform distribution from 0.0 to 1.0
};

} /* namespace networkVMC */

#endif /* SRC_UTILITIES_RNGWRAPPER_HPP_ */
