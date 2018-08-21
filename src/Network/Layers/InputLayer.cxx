/*
 * InputLayer.cxx
 * Created on 06.3.2018
 * Author: Ke Liao 
 */

#include <vector>
#include <Eigen/Dense>
#include "../../HilbertSpace/Determinant.hpp"
#include "InputLayer.hpp"

namespace networkVMC{

template <typename F, typename coeffType>
InputLayer<F, coeffType>::InputLayer(std::vector<Eigen::VectorXd> const &inputs_, int size_):
    Layer(inputs_, "Linear"), numNrn(size_){
  //for the InputLayer, the inputs are 0 vector.
  activations.resize(1,Eigen::VectorXd::Zero(numNrn));
}

template <typename F, typename coeffType>
InputLayer<F, coeffType>::~InputLayer(){
}

}
