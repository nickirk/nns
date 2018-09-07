/*
 * RBM.hpp
 *
 *  Created on: Apr 29, 2018
 *      Author: ghanem and Liao
 */
#ifndef RBM_DEFINED
#define RBM_DEFINED

#include <Eigen/Dense>
#include "Parametrization.hpp"
#include "../utilities/StateForward.hpp"
#include "../utilities/TypeDefine.hpp"

namespace networkVMC{
  //Resttricted Boltzmann Machine
  template <typename F=std::complex<double>, typename coeffType=std::complex<double>>
  class RBM: public ClonableParametrization<F, coeffType, RBM<F, coeffType>>{
  public:
    using T=Eigen::Matrix<F, Eigen::Dynamic, 1>;
    RBM(int sizeInput_, int sizeHidden_);

    coeffType getCoeff(detType const &det) const;
    T getDeriv(detType const &det) const;
    T const& pars() const;

    // copying has to be defined explicitly, because it requires re-initializing the maps
    RBM(RBM<F, coeffType> const &source);
    RBM<F, coeffType> & operator=(RBM<F, coeffType> const &);
    // but moving works
    RBM(RBM<F, coeffType> &&) = default;
    RBM<F, coeffType> & operator=(RBM<F, coeffType> &&) = default;
    ~RBM();

  private:
    int sizeHidden;
    int sizeInput;
    int numPars;
    int a_offset;
    int b_offset;
    int w_offset;
    T pars_vec;
    Eigen::Map<T> a;
    Eigen::Map<T> b;
    Eigen::Map<Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>> w;
    };
}

#endif //RBM_DEFINED
