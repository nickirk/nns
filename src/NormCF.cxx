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
	calced = false;
	double buf{0.0};
	for(int i=0;i<input.size();++i){
		coeffType tmp = (input.getCoeff(i) - psi.getCoeff(i));
		buf += tmp.sum();
	}
	return 2*buf;
}

double NormCF::calc(State const &input) const{
	double buf{0.0};
	for(int i=0;i<input.size();++i){
		coeffType tmp = input.getCoeff(i) - psi.getCoeff(i);
		coeffType ctmp = tmp;
		ctmp[1] = -ctmp[1];
		buf += ctmp*tmp;
	}
	return buf;
}


