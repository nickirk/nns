/*
 * calcNabla.hpp
 *
 *  Created on: Aug 23, 2018
 *      Author: guther
 */

#ifndef SRC_NETWORK_CALCNABLA_HPP_
#define SRC_NETWORK_CALCNABLA_HPP_

#include "../utilities/TypeDefine.hpp"
#include "../utilities/MatrixTypeDefine.hpp"
#include "Parametrization.hpp"
#include "../utilities/State.hpp"

#include <Eigen/Dense>

namespace networkVMC{

template<typename D>
paraVector calcNabla(State const &inputState,
		paraVector const &dEdC, paraType const &energy, D &&deriv, SamplerType sT, int numPars){

	int numDets = inputState.size();
	int spaceSize = inputState.totalSize();
	paraVector dEdW= paraVector::Zero(numPars);
	Eigen::Matrix<paraType, Eigen::Dynamic, Eigen::Dynamic>
	dCdW = Eigen::Matrix<paraType, Eigen::Dynamic, Eigen::Dynamic>::Zero(numPars, spaceSize);
#pragma omp parallel
	{
	// fill up the matrix of dCdW, like in EnergyEsMarkov.cxx
	// reserve space and in the end use matrix*vector instead of
	// a summation
	paraVector dCtdW;
#pragma omp for
	for (int i=0; i < numDets; ++i){
		//need private dCtdW
		// do the mapping inside for loop, private
		//update vector dCidWk
		dCtdW = deriv(inputState.det(i));
		// multiplication should be done by matrix vector product
		// fill up the dCdW matrix
		dCdW.col(i) << (dCtdW.conjugate());
		if(sT==Markov) dEdW -= energy * dCtdW.conjugate()/ inputState.getTotalWeights();
		//dedc[i] = 1;
		std::vector<detType> coupledDets = inputState.coupledDets(i);
		std::vector<coeffType> coupledCoeffs = inputState.coupledCoeffs(i);
		size_t coupledSize = inputState.coupledDets(i).size();
		size_t pos = inputState.locate(i);
		for (size_t j(0); j < coupledSize; ++j){
			// fill up the dCdW matrix with coupled dets contribution
		  //dCtdW = deriv(coupledDets(j]);
			dCdW.col(numDets+pos+j) << (dCtdW.conjugate());
			//dEdWTmp +=  dCtdW * dEdC[pos];
		}
	}
	}
	dEdW += (dCdW * dEdC);//.conjugate();
	return dEdW;
}

} // end namespace

#endif /* SRC_NETWORK_CALCNABLA_HPP_ */
