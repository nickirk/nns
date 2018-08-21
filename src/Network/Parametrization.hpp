/*
 * Parametrization.hpp
 *
 *  Created on: Apr 30, 2018
 *      Author: guther
 */

#ifndef SRC_NETWORK_PARAMETRIZATION_HPP_
#define SRC_NETWORK_PARAMETRIZATION_HPP_

#include "../utilities/TypeDefine.hpp"
#include <fstream>
#include <iterator>
#include <string>
#include <utility>
#include <Eigen/Dense>
#include "../utilities/State.hpp"
#include "../utilities/Errors.hpp"

namespace networkVMC {

// Forward declaration for more efficient compilation
//class State<double>;
//class State<std::complex<double>>;

// This is the base class for any parametrization of the wave function
// It is an abstract parametrization that can be updated to be optimized
// with respect to a given cost function
template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
class Parametrization {
  public:
	using T=Eigen::Matrix<F, Eigen::Dynamic ,1>;
    Parametrization(){};
    // default the big five - required for virtual destructor
    virtual ~Parametrization() = default;

    // It needs to be able to return coefficients somehow
    virtual coeffType getCoeff(detType const &det) const=0; // Can throw an invalidDeterminantError
    // base method for returning the parameters (as a single vector)
    virtual T const& pars() const = 0;
    // We also need write access to the parameters (this is the non-const variant)
    virtual T& pars(){return const_cast<T&>
      (static_cast<Parametrization<F, coeffType> const&>(*this).pars());};
    virtual void writeParsToFile(std::string file){
      T tmpPars = pars();
      std::vector<F> parsStd(tmpPars.data(), tmpPars.data()+tmpPars.size());
      std::ofstream fout(file);
      fout.precision(10);
      std::copy(parsStd.begin(), parsStd.end(),
      std::ostream_iterator<F>(fout, "\n"));
      fout.close();
    };
    virtual void readParsFromFile(std::string file){
      std::ifstream input(file);
      if (!input){
         throw FileNotFound(file);
         abort();
      }

      T &innerPars = pars();
      std::vector<F> buffer{
        std::istream_iterator<F>(input),
        std::istream_iterator<F>() };
      assert (buffer.size() == pars().size());
      innerPars = Eigen::Map<T> (buffer.data(),buffer.size());
      //use the pars() with write access
    };
    virtual int getNumPars(){ return pars().size();};
//   Obtain the inner derivative dX/dPars with given dX/dC (C are coefficients)
    /*
    virtual T calcNablaPars(
  		  State const &input,
  		  T const &outerDerivative) = 0; // can throw a SizeMismatchError if input and outerDerivative
    	  	  	  	  	  	  	  	  	  	  	  	 // have different size
    */
    virtual T getMarkovDeriv(detType const &det) const{
  	  return getDeriv(det)/getCoeff(det);
    };
    virtual T getDeriv(detType const &det) const{
  	  return T();
    };
//   The following features are experimental and not essential to the interface
    // Some other derivative
  // The following function is used when use ListGen or fullSampler
  template <typename F, typename coeffType>
  T calcNablaParsConnected(
    State<coeffType> const &inputState,
    T const &dEdC
   ){
      int numDets = inputState.size();
      // spaceSize = size of sampled dets and their coupled ones
      int spaceSize = inputState.totalSize();
      int numPars=getNumPars();
      T dEdW= T::Zero(numPars);
      // T dEdWTmp= T::Zero(numPars);
      Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> dCdW = Eigen::MatrixXcd::Zero(numPars, spaceSize);
      std::vector<std::complex<double>> dedc=dEdC;
        #pragma omp for
        // fill up the matrix of dCdW, like in EnergyEsMarkov.cxx
        // reserve space and in the end use matrix*vector instead of
        // a summation
        for (int i=0; i < numDets; ++i){
          //need private dCtdW
          T dCtdW= T::Zero(numPars);
          // do the mapping inside for loop, private
          //update vector dCidWk
          dCtdW = getDeriv(inputState.det(i));
          //update vector dCidWk
          // multiplication should be done by matrix vector product
          // fill up the dCdW matrix
          dCdW.col(i) << (dCtdW.conjugate());
          //dedc[i] = 1;
          std::vector<detType> coupledDets = inputState.coupledDets(i);
          std::vector<coeffType> coupledCoeffs = inputState.coupledCoeffs(i);
          size_t coupledSize = inputState.coupledDets(i).size();
          size_t pos = inputState.locate(i);
          for (size_t j(0); j < coupledSize; ++j){
            //update dCjdW
            //dCtdW = getDeriv(coupledDets[j]);
            // fill up the dCdW matrix with coupled dets contribution
            dCdW.col(numDets+pos+j) << (dCtdW.conjugate());
            //dEdWTmp +=  dCtdW * dEdC[pos];
          }
        }
      // map std::vector of dEdC to Eigen Vector
      T dEdCEigen=Eigen::Map<T>(dedc.data(),spaceSize);
      // make it parallel. TODO
      dEdW = (dCdW * dEdCEigen);//.conjugate();
      return dEdW;
    }

    virtual T calcNablaParsMarkovConnected(State<coeffType> const &inputState,
                  T const& dEdC, F const& energy){
      int numDets = inputState.size();
      int spaceSize = inputState.totalSize();
      int numPars=getNumPars();
      T dEdW= T::Zero(numPars);
      Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>
      dCdW = Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>::Zero(numPars, spaceSize);
      #pragma omp for
      // fill up the matrix of dCdW, like in EnergyEsMarkov.cxx
      // reserve space and in the end use matrix*vector instead of
      // a summation
      for (int i=0; i < numDets; ++i){
        //need private dCtdW
        T dCtdW= T::Zero(numPars);
        // do the mapping inside for loop, private
        //update vector dCidWk
        dCtdW = getMarkovDeriv(inputState.det(i));
        // multiplication should be done by matrix vector product
        // fill up the dCdW matrix
        dCdW.col(i) << (dCtdW.conjugate());
        dEdW -= energy * dCtdW.conjugate()/numDets;
        //dedc[i] = 1;
        std::vector<detType> coupledDets = inputState.coupledDets(i);
        std::vector<coeffType> coupledCoeffs = inputState.coupledCoeffs(i);
        size_t coupledSize = inputState.coupledDets(i).size();
        size_t pos = inputState.locate(i);
        for (size_t j(0); j < coupledSize; ++j){
          // fill up the dCdW matrix with coupled dets contribution
          dCdW.col(numDets+pos+j) << (dCtdW.conjugate());
          //dEdWTmp +=  dCtdW * dEdC[pos];
        }
      }
      // map std::vector of dEdC to Eigen Vector
      //T dEdCEigen=Eigen::Map<T>(dedc.data(),spaceSize);
      // make it parallel. TODO
      //dEdW += (dCdW * dEdCEigen);//.conjugate();
      dEdW += (dCdW * dEdC);//.conjugate();
      return dEdW;
    }

    // stochastic reconfiguration derivative
    /*
    virtual Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> calcdCdwSR(
      State const &outputState
    ){return Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>();};
    */
    // virtual construction
    virtual Parametrization<F, coeffType>* clone() const = 0;
    virtual Parametrization<F, coeffType>* move_clone() = 0;

};

// implement the polymorphic copy for the derived classes
template <typename T, typename coeffType, typename Base>
class ClonableParametrization: public Parametrization<T, coeffType>{
  public:
	// inherit the constructor
	using Parametrization<T, coeffType>::Parametrization;
	virtual ~ClonableParametrization() = default;

	// copy-clone (does not change *this)
	virtual Parametrization<T, coeffType>* clone() const{
		return new Base{static_cast<Base const&>(*this)};
	}

	// move-clone (moves the data to the return pointer)
	virtual Parametrization<T, coeffType>* move_clone(){
		return new Base(std::move(static_cast<Base &>(*this) ) );
	}
};


} /* namespace networkVMC */

#endif /* SRC_NETWORK_PARAMETRIZATION_HPP_ */
