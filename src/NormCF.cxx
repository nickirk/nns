/*
 * NormCF.cxx
 *
 *  Created on: Oct 26, 2017
 *      Author: guther
 */

#include <math.h>
#include "NormCF.hpp"
#include "CoeffType.hpp"

std::vector<coeffType > NormCF::nabla(State const &input) const{
	std::vector<coeffType > cfBuf(input.size());
	for(size_t i=0;i<input.size();++i){
		coeffType buf;
		buf = Eigen::VectorXd::Zero(2);
		for(size_t i=0;i<input.size();++i){
			buf += (input.getCoeff(i) - psi.getCoeff(i));
		}
		cfBuf[i] = buf;
	}
	return cfBuf;
}

double NormCF::calc(State const &input) const{
	double buf{0.0};
	for(size_t i=0;i<input.size();++i){
		coeffType tmp = input.getCoeff(i) - psi.getCoeff(i);
		buf += pow(tmp[0],2) + pow(tmp[1],2);
	}
	return buf;
}


