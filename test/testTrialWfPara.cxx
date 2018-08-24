/*
 * testTrialWfPara.cxx
 *
 *  Created on: Aug 16, 2018
 *      Author: guther
 */

#include <iostream>
#include <assert.h>
#include <complex>
#include <Eigen/Dense>
#include "../src/NNWLib.hpp"
#include "defaultSystem.hpp"

using namespace networkVMC;

int main(){
	int numSites = 4;
	int nHidden = 8;
	auto twfBasis = generateDefaultBasis(numSites);
	// do a simple test of the TrialWfPara with
	// an RBM as a base para
	RBM<> network(numSites*2,nHidden);
	auto HF = twfBasis.getDetByIndex(0);

	// and a linear parametrization as a trial WF
	auto trial = DirectParametrization<std::complex<double>>(twfBasis);

	TrialWfPara<std::complex<double>> twfPara(network,trial);

	// first compare the coefficients of the base*trial with that of the TrialWfPara
	assert(std::abs(twfPara.getCoeff(HF) - twfPara.getTrialCoeff(HF)*twfPara.getBaseCoeff(HF)) < epsilon);

  std::cout << "twfPara.getBaseCoeff(HF)= " << twfPara.getBaseCoeff(HF) << std::endl;
  std::cout << "twfPara.getTrialCoeff(HF)= " << twfPara.getTrialCoeff(HF) << std::endl;
  std::cout << "twfPara.getCoeff(HF)= " << twfPara.getCoeff(HF) << std::endl;
  std::cout << "twfPara.getBaseCoeff(HF)*twfPara.getTrialCoeff(HF)= " << 
    twfPara.getBaseCoeff(HF)*twfPara.getTrialCoeff(HF) << std::endl;
	// and now check the derivative (w.r. to 1 coefficient)

	// define a dummy state with 1 det and 3 times the same det as coupled (all with coeff. of 1)
	State<> singleDet(1);
	singleDet.coeff(0) = 1.0;
	singleDet.det(0) = HF;
	singleDet.coupledDets(0) = std::vector<detType>(0,HF);
	singleDet.coupledCoeffs(0) = std::vector<cType>(0,1.0);

	Eigen::VectorXcd outerDerivative(1);
  outerDerivative << (cType(1,1.0));

	auto nabla = twfPara.calcNablaParsConnected(singleDet,outerDerivative);

	// the base parametrization
	Parametrization<cType, cType> const& baseRef = twfPara.base();
	// the number of parameters does not change
	assert(nabla.size() == baseRef.pars().size());

	// simulate polymorphic access
	Parametrization<cType, cType> &twfRef = twfPara;
	// now test manipulation of the trial wf pars
	std::cout << "Parameters checksum: "<< twfRef.pars().squaredNorm() << std::endl;
	twfRef.pars() *= std::complex<double>(2.0,0.0);
	std::cout << "Parameters checksum: "<< twfRef.pars().squaredNorm() << std::endl;
	std::cout << "Coeff now: " << twfRef.getCoeff(HF) << std::endl;
	std::cout << "Base checksum: " << baseRef.pars().squaredNorm() << std::endl;
	std::cout << "Base coeff: " << baseRef.getCoeff(HF) << std::endl;
}


