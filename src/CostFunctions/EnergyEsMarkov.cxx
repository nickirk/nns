/*
 * EnergyEsMarkov.cxx
 * based on EnergyCF.cxx
 *  Created on: Jan 29, 2018
 *      Author: Liao
 */

#include "EnergyEsMarkov.hpp"

#include <vector>
#include <complex>
#include <iostream>

#include "../utilities/State.hpp"
#include "../Hamiltonian/Hamiltonian.hpp"

namespace networkVMC{

coeffType EnergyEsMarkov::evaluate(State const &input) const{
  coeffType energyVal{0.0};
  normalizerCoeff=0.0;
  int numDets = input.size();
  //#pragma omp parallel for reduction(+:energyVal)
  for (int i=0; i < numDets; ++i){
    double Hij(0.);
    coeffType c_i=input.coeff(i);
    std::vector<coeffType> coupledC_j = input.coupledCoeffs(i);
    std::vector<detType> coupledDets = input.coupledDets(i);
    std::vector<double> coupledWeights = input.coupledWeights(i);
    Hij = H(input.det(i), input.det(i));
    energyVal += c_i * Hij/c_i;
    //std::cout << "ci=" << c_i << std::endl;
    for (size_t j=0; j < coupledC_j.size(); ++j){
        // don't forget to unbias using the Pgen. TODO
      Hij = H(input.det(i), coupledDets[j]);
      energyVal += coupledC_j[j]* Hij / c_i/
          coupledC_j.size()/coupledWeights[j];
    }
  }
  energyVal /= numDets;
  return energyVal;
}

nablaType EnergyEsMarkov::nabla(State const &input) const{
  coeffType energyM = evaluate(input);
  energy = std::real(energyM);
  int numDets = input.size();
  // spaceSize = size of sampled dets and their coupled ones
  int spaceSize = input.totalSize();
  // For each n, it has the same number of connected Dets as others.!!!
  nablaType dEdC(spaceSize, coeffType(0.,0.)) ;
  //not thread safe
  //assume we know the whole space size, reserve space
  #pragma omp parallel for
  for (int i=0; i < numDets; ++i){
    coeffType dEdCtmp;
    coeffType c_i = input.coeff(i);
//  put all the weighting step here instead of inside of RBM
    dEdCtmp = (H(input.det(i), input.det(i)) - energyM)/c_i;
    // add weights
    dEdC[i] = dEdCtmp / numDets;
    std::vector<coeffType> coupledC_j = input.coupledCoeffs(i);
    std::vector<detType> coupledDets = input.coupledDets(i);
    std::vector<double> coupledWeights = input.coupledWeights(i);
    int tid = omp_get_thread_num();
    std::cout << "numProc=" << tid << " EnergyEsMarkov.cxx: c_i=" << c_i << std::endl;
    int coupledSize = input.coupledDets(i).size();
   int pos=input.locate(i);
   for (size_t j=0; j < coupledSize; ++j){
     // don't forget to unbias using the Pgen. TODO
     // in the excitgen, the weight should be updated with pgen
     dEdCtmp = H(input.det(i),coupledDets[j])/c_i/ coupledC_j.size()/coupledWeights[j];
     //std::cout << "EnergyEsMarkov.cxx: coupledC_j=" << coupledC_j[j] << std::endl;
     //std::cout << "EnergyEsMarkov.cxx: coupledC_j size=" << coupledC_j.size() << std::endl;
     //std::cout << "EnergyEsMarkov.cxx: coupledWeights j=" << coupledWeights[j] << std::endl;
     // unbias with numCoupledDets and Pgen
     //dEdC[i+j+1]=dEdCtmp / coupledDets.size() / numDets;
      dEdC[numDets+pos+j]=dEdCtmp / numDets;
    }
  }
  // dEdC is stored as folllowing:
  // |c_0|c_1|...|c_{numDets-1}|c_0^{0}|c_0^{1}|...|c_0^{coupledSize-1}|c_1^{0}
  // |c..|c_{numDets-1}^{0}|...|c_{numDets-1}^{coupledSize-1}|
  return dEdC;
}

}

