/*
 * NormCF.cxx
 *
 *  Created on: Oct 26, 2017
 *      Author: guther
 */

#include "NormCF.hpp"

#include <math.h>

namespace networkVMC{

nablaType NormCF::nabla(State const &input) const{
	nablaType cfBuf(input.size());
	for(std::size_t i=0;i<input.size();++i){
		coeffType buf0;
		Eigen::VectorXd buf = Eigen::VectorXd::Zero(2);
		for(size_t j=0;j<input.size();++j){
			buf[0] += std::real(input.coeff(j) - psi.coeff(j));
			buf[1] += std::imag(input.coeff(j) - psi.coeff(j));
		}
		buf0.real(buf[0]);
		buf0.imag(buf[1]);
		cfBuf[i] = 2*buf0;
	}
	return cfBuf;
}

double NormCF::calc(State const &input) const{
	double buf{0.0};
	for(std::size_t i=0;i<input.size();++i){
		coeffType tmp = input.coeff(i) - psi.coeff(i);
		buf += std::real(std::conj(tmp) * tmp);
	}
	return buf;
}

}
