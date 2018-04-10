/*
 * OutputLayer.cxx
 * Created on 10.04.2018
 * Author: Ke Liao 
 */


#include <vector>
#include <Eigen/Dense>
#include "OutputLayer.hpp"


OutputLayer(std::vector<Eigen::VectorXd> const &inputs_, 
           double &(actFunc_)(double), int size_):
  DenseLayer(inputs_, &(actFunc_)(double), size_);


void DenseLayer::backProp(
    std::vector<Eigen::VectorXd> prevDelta, 
    weightType &prevWeights;
    ){
    deltas[0] = prevWeights[0][0].transpose() * prevDelta; 
    deltas[0] = deltas[0].array()* z[0].unaryExpr(&actFuncPrime).array();
    nablaBiases[0] = deltas[0];
    //get a weight matrix. \partial C/\partial w^{l}_{jk} = a^{l-1}_k \delta_j^l 
    //the layer here refers to the lth layer of Biases and weights, so for
    //activation layer refers to the l-1th layer.
    for (size_t i(0); i < inputs.size(); i++)
      nablaWeights[0][i] = delta[0] * inputs[i].transpose();
}
