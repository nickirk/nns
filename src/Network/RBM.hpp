/*
 * RBM.hpp
 *
 *  Created on: Apr 29, 2018
 *      Author: ghanem and Liao
 */
#ifndef RBM_DEFINED
#define RBM_DEFINED

#include <Eigen/Dense>
#include "../CostFunctions/CostFunction.hpp"
#include "Parametrization.hpp"
#include "../utilities/State.hpp"
#include "../utilities/TypeDefine.hpp"
#include "../math/MathFunctions.hpp"

namespace networkVMC{
  //Resttricted Boltzmann Machine
  template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
  class RBM: public ClonableParametrization<F, coeffType, RBM<F, coeffType>>{
  public:
    using T=Eigen::Matrix<F, Eigen::Dynamic, 1>;
    RBM(int sizeInput_, int sizeHidden_);

    coeffType getCoeff(detType const &det) const;
    T getDeriv(detType const &det) const;
    T getMarkovDeriv(detType const &det) const;
    T const& pars() const;
    //T calcNablaPars(State<coeffType> const &input, T const &outerDerivative);
    //Eigen::Matrix<F, Dynamic, Dynamic> calcdCdwSR(State const &outputState);
    //T calcNablaParsConnected(State<coeffType> const &inputState, T const& dEdC);
    //T calcNablaParsMarkovConnected(State const &inputState, T const& dEdC, F const& energy);

  private:
    int sizeHidden;
    int numPars;
    int sizeInput;
    int a_offset;
    int b_offset;
    int w_offset;
    T pars_vec;
    Eigen::Map<T> a;
    Eigen::Map<T> b;
    Eigen::Map<T> w;
    };
}

#endif //RBM_DEFINED
