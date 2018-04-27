/*
 * NormCF.cxx
 *
 *  Created on: Oct 26, 2017
 *      Author: guther
 */

#include "NormCF.hpp"

#include <math.h>

#include "../utilities/TypeDefine.hpp"

namespace networkVMC{

std::vector<Eigen::VectorXd > NormCF::nabla(State const &input) const{
	std::vector<Eigen::VectorXd > cfBuf(input.size());
	for(size_t i=0;i<input.size();++i){
		Eigen::VectorXd buf;
		buf = Eigen::VectorXd::Zero(2);
		for(size_t j=0;j<input.size();++j){
			buf[0] += std::real(input.coeff(j) - psi.coeff(j));
			buf[1] += std::imag(input.coeff(j) - psi.coeff(j));
		}
		cfBuf[i] = 2*buf;
	}
	return cfBuf;
}

double NormCF::calc(State const &input) const{
	double buf{0.0};
	for(size_t i=0;i<input.size();++i){
		coeffType tmp = input.coeff(i) - psi.coeff(i);
		buf += std::real(std::conj(tmp) * tmp);
	}
	return buf;
}

}

