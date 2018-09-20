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
#include "../utilities/TypeDefine.hpp"

namespace networkVMC{

class State;

  //Resttricted Boltzmann Machine
  class RBM: public ClonableParametrization<RBM>{
  public:
    RBM(int sizeInput_, int sizeHidden_);

    coeffType getCoeff(detType const &det) const;
    paraVector getDeriv(detType const &det) const;
    paraVector const& pars() const;

    // copying has to be defined explicitly, because it requires re-initializing the maps
    RBM(RBM const &source);
    RBM & operator=(RBM const &);
    // but moving works
    RBM(RBM &&) = default;
    RBM & operator=(RBM &&) = default;
    ~RBM();

  private:
    int sizeHidden;
    int sizeInput;
    int numPars;
    int a_offset;
    int b_offset;
    int w_offset;
    paraVector pars_vec;
    Eigen::Map<paraVector> a;
    Eigen::Map<paraVector> b;
    Eigen::Map<Eigen::Matrix<paraType, Eigen::Dynamic, Eigen::Dynamic>> w;
    };
}

#endif //RBM_DEFINED
