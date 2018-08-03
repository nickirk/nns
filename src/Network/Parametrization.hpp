/*
 * Parametrization.hpp
 *
 *  Created on: Apr 30, 2018
 *      Author: guther
 */

#ifndef SRC_NETWORK_PARAMETRIZATION_HPP_
#define SRC_NETWORK_PARAMETRIZATION_HPP_

#include "../utilities/TypeDefine.hpp"

namespace networkVMC {

// Forward declaration for more efficient compilation
class State;

// This is the base class for any parametrization of the wave function
// It is an abstract parametrization that can be updated to be optimized
// with respect to a given cost function
template<typename T=VecType>
class Parametrization {
public:
  Parametrization(){};
  virtual ~Parametrization(){};
 // It needs to be able to return coefficients somehow
  virtual coeffType getCoeff(detType const &det) const=0; // Can throw an invalidDeterminantError
  // base method for returning the parameters (as a single vector)
  virtual T const& pars() const = 0;
  // We also need write access to the parameters (this is the non-const variant)
  virtual T& pars(){return const_cast<T&>
    (static_cast<Parametrization const&>(*this).pars());};
  virtual int getNumPars(){ return pars().size();};
// Obtain the inner derivative dX/dPars with given dX/dC (C are coefficients)
  virtual T calcNablaPars(
		  State const &input,
		  nablaType const &outerDerivative) = 0;

  virtual Eigen::VectorXcd getMarkovDeriv(detType const &det) const{
	  return getDeriv(det)/getCoeff(det);
  };
  virtual Eigen::VectorXcd getDeriv(detType const &det) const{
	  return Eigen::VectorXcd();
  };
// The following features are experimental and not essential to the interface
  // Some other derivative
  virtual T calcNablaParsConnected(State const &inputState, nablaType const& dEdC)
  	  {return calcNablaPars(inputState,dEdC);}

  virtual VecCType calcNablaParsMarkovConnected(State const &inputState, 
                nablaType const& dEdC, coeffType const& energy){
    int numDets = inputState.size();
    int spaceSize = inputState.totalSize();
    int numPars=getNumPars();
    Eigen::VectorXcd dEdW= Eigen::VectorXcd::Zero(numPars);
    Eigen::MatrixXcd dCdW = Eigen::MatrixXcd::Zero(numPars, spaceSize);
    std::vector<std::complex<double>> dedc=dEdC;
    #pragma omp for
    // fill up the matrix of dCdW, like in EnergyEsMarkov.cxx
    // reserve space and in the end use matrix*vector instead of
    // a summation
    for (int i=0; i < numDets; ++i){
      //need private dCtdW
      Eigen::VectorXcd dCtdW= Eigen::VectorXcd::Zero(numPars);
      // do the mapping inside for loop, private
      //update vector dCidWk
      dCtdW = getMarkovDeriv(inputState.det(i));
      // multiplication should be done by matrix vector product
      // fill up the dCdW matrix
      dCdW.col(i) << (dCtdW.conjugate());
      dEdW -= energy * dCtdW.conjugate()/numDets; 
      //dedc[i] = 1;
      std::vector<detType> coupledDets = inputState.coupledDets(i);
      std::vector<coeffType > coupledCoeffs = inputState.coupledCoeffs(i);
      size_t coupledSize = inputState.coupledDets(i).size();
      size_t pos = inputState.locate(i);
      for (size_t j(0); j < coupledSize; ++j){
        // fill up the dCdW matrix with coupled dets contribution
        dCdW.col(numDets+pos+j) << (dCtdW.conjugate());
        //dEdWTmp +=  dCtdW * dEdC[pos];
      }
    }
    // map std::vector of dEdC to Eigen Vector
    Eigen::VectorXcd dEdCEigen=Eigen::Map<Eigen::VectorXcd>(dedc.data(),spaceSize);
    // make it parallel. TODO
    dEdW += (dCdW * dEdCEigen);//.conjugate();
    return dEdW;
  }

  // stochastic reconfiguration derivative
  virtual Eigen::MatrixXcd calcdCdwSR(
    State const &outputState
  ){return Eigen::MatrixXd();};

};


} /* namespace networkVMC */

#endif /* SRC_NETWORK_PARAMETRIZATION_HPP_ */
