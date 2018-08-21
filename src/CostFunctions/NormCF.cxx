/*
 * NormCF.cxx
 *
 *  Created on: Oct 26, 2017
 *      Author: guther
 */

#include "NormCF.hpp"

#include <math.h>

namespace networkVMC{
template <typename F, typename coeffType>
NormCF<F, coeffType>::T NormCF<F, coeffType>::nabla(State<coeffType> const &input) const{
	T cfBuf(input.size());
	for(std::size_t i=0;i<input.size();++i){
		coeffType buf0;
		Eigen::VectorXd buf = Eigen::VectorXd::Zero(2);
		//buf[0] = std::real(input.coeff(i) - psi.coeff(i));
		//buf[1] = std::imag(input.coeff(i) - psi.coeff(i));
    buf0 = std::conj(input.coeff(i)-psi.coeff(i));
		//buf0.real(buf[0]);
		//buf0.imag(0.);
		cfBuf[i] = 2.0*buf0;
	}
	return cfBuf;
}

template <typename F, typename coeffType>
coeffType NormCF<F, coeffType>::calc(State<coeffType> const &input) const{
	coeffType buf=0.;
	for(std::size_t i=0;i<input.size();++i){
		coeffType tmp = input.coeff(i) - psi.coeff(i);
		buf += std::real(std::conj(tmp) * tmp);
	}
	return buf;
}
template class NormCF<double, double>;
template class NormCF<std::complex<double>, std::complex<double>>;

}

