/*
 * AcceleratedGradientDescent.cxx
 *
 *  Created on: May 2, 2018
 *      Author: guther
 */

#include "AcceleratedGradientDescent.hpp"
#include "../utilities/Errors.hpp"

namespace networkVMC {

template <typename F, typename coeffType>
AcceleratedGradientDescent<F, coeffType>::AcceleratedGradientDescent(double learningRate_):
  Solver<F, coeffType>(learningRate_),numPars(0),lambdaS1(0.0),lambdaS(0.0),gammaS(0.0),
  gammaS1(0.0),uninitialized(true){
    yS = T::Zero(numPars);
    yS1 = T::Zero(numPars);
    Egz2 = T::Zero(numPars);
}

//---------------------------------------------------------------------------------------------------//

template <typename F, typename coeffType>
AcceleratedGradientDescent<F, coeffType>::~AcceleratedGradientDescent() {
}

//---------------------------------------------------------------------------------------------------//

// This is ported from Nnw.cxx
template <typename F, typename coeffType>
void AcceleratedGradientDescent<F, coeffType>::update(T &w, T const &force,
		State<coeffType> const &input, SamplerType const & samplerType){
  if(uninitialized){
	numPars = w.size();
    yS = T::Zero(numPars);
    yS1 = T::Zero(numPars);
    Egz2 = T::Zero(numPars);
  }
  // we can only proceed, if the same number of parameters is given in each update
  if(w.size()!=numPars) throw errors::SizeMismatchError(w.size(),numPars);

  // Do the update
  lambdaS1 = (1+std::sqrt(1+4*lambdaS*lambdaS))/2.;
  gammaS = (1-lambdaS)/(lambdaS1);//*std::exp(-1./100*iteration);
  double rho=0.9;//*std::exp(-1./100*iteration);
  Egz2 = rho*Egz2.matrix()+(1-rho)*force.array().square().matrix();
  Egz2 += T::Ones(numPars)*1e-4;
  T RMS = (Egz2).array().sqrt();
  T tau = Solver<F, coeffType>::learningRate * RMS.array().inverse();
  //yS1 = NNP - learningRate * generlisedForce;
  yS1 = (w.array() -  tau.array() * force.array()).matrix();

  // if this is the first update, note this
  if (uninitialized){
	  yS=yS1;
	  uninitialized = false;
  }
  // set the output
  w = (1-gammaS)*yS1 + gammaS*yS;
  lambdaS = lambdaS1;
  yS = yS1;
}
//instantiate class
template class AcceleratedGradientDescent<double, double>;

template class AcceleratedGradientDescent<std::complex<double>, std::complex<double>>;
} /* namespace networkVMC */
