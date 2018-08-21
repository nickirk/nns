/*
 * ADAM.cxx
 *
 *  Created on: May 2, 2018
 *      Author: guther
 */

#include "ADAM.hpp"

namespace networkVMC {
template <typename F, typename coeffType>
ADAM<F, coeffType>::ADAM(double learningRate_):Solver<F, coeffType>(learningRate_),numPars(0),
		iteration(0),beta1(0.9),beta2(0.999) {
	m  = T::Zero(numPars);
	m1  = T::Zero(numPars);
	v  = T::Zero(numPars);
	v1  = T::Zero(numPars);
}

template <typename F, typename coeffType>
ADAM<F, coeffType>::~ADAM() {
}

template <typename F, typename coeffType>
void ADAM<F, coeffType>::update(T &w, T const &force, State<coeffType> const &input, SamplerType const
    &samplerType){
	// in the first iteration, we learn about the number of parameters
	if(iteration==0){
		numPars = w.size();
		m  = T::Zero(numPars);
		m1  = T::Zero(numPars);
		v  = T::Zero(numPars);
		v1  = T::Zero(numPars);
	}
	iteration+=1;
    m1 = beta1*m + (1-beta1)*force;
    v1 = beta2*v + (1-beta2)*force.array().square().matrix();
    m = m1;
    v = v1;
    m1 = m1/(1-std::pow(beta1,iteration));
    v1 = v1/(1-std::pow(beta2,iteration));
    w -= Solver<F, coeffType>::learningRate * (m1.array()/(v1.array().sqrt()+1e-8)).matrix();
}
//instatiate class
template class ADAM<std::complex<double>, std::complex<double>>;
;
template class ADAM<double, double>;
} /* namespace networkVMC */
