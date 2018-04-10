/*
 * NormCF.cxx
 *
 *  Created on: Oct 26, 2017
 *      Author: guther
 */

#include <math.h>
#include "NormCF.hpp"
#include "TypeDefine.hpp"

std::vector<Eigen::VectorXd > NormCF::nabla(std::vector<State> const &input) const{
	std::vector<Eigen::VectorXd > cfBuf(input.size());
	for(size_t i=0;i<input.size();++i){
		Eigen::VectorXd buf;
		buf = Eigen::VectorXd::Zero(2);
		for(size_t j=0;j<input.size();++j){
			buf[0] += std::real(input[j].coeff - psi[j].coeff);
			buf[1] += std::imag(input[j].coeff - psi[j].coeff);
		}
		cfBuf[i] = 2*buf;
	}
	return cfBuf;
}

double NormCF::calc(std::vector<State> const &input) const{
	double buf{0.0};
	for(size_t i=0;i<input.size();++i){
		coeffType tmp = input[i].coeff - psi[i].coeff;
		buf += std::real(std::conj(tmp) * tmp);
	}
	return buf;
}


