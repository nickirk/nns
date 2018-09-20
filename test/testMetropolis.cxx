/*
 * testMetropolis.cxx
 *
 *  Created on: Apr 30, 2018
 *      Author: guther
 */

#include "../src/NNWLib.hpp"
#include "defaultSystem.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <unordered_map> //for unordered_map
//TODO : Write a statistical test that checks if the correct probability
//       density is generated

using namespace networkVMC;

double testL2Norm(std::vector<int> directSample,
    std::vector<int> markovSample, int nSweeps){
  //make sure the two vectors are of the same size
  assert (directSample.size()==markovSample.size());
  double Z(0.);
  //write data to file for visualisation
  std::ofstream myfile;
  myfile.open("Dist.dat");
  for(int i(0); i < directSample.size(); ++i) {
    //std::cout << "Frequency direct and markov=" << directSample[i] << " " << markovSample[i]
    //  << std::endl;
    myfile << i << "  " << directSample[i]/double(nSweeps) << "  " << 
      markovSample[i]/double(nSweeps) << std::endl;
    Z += double(std::pow(directSample[i]-markovSample[i],2)-directSample[i]-markovSample[i]);
  }
  myfile.close();
  return std::sqrt(Z)/nSweeps;
}

double testL1Norm(std::vector<int> directSample,
    std::vector<int> markovSample, int nSweeps){
  //make sure the two vectors are of the same size
  assert (directSample.size()==markovSample.size());
  double Z(0.);
  //write data to file for visualisation
  for(int i(0); i < directSample.size(); ++i) {
    //std::cout << "Frequency direct and markov=" << directSample[i] << " " << markovSample[i]
    //  << std::endl;
    if ((directSample[i]+markovSample[i])==0) continue;
    Z += double(std::pow(directSample[i]-markovSample[i],2)-directSample[i]-markovSample[i])
         / double(directSample[i]+markovSample[i]);
  }
  Z /= std::sqrt(nSweeps);
  //std::REQUIRE(Z <= 5.);
  return Z;
}
std::vector<int> genDirectSample(Parametrization
    const &para, Basis const &basis, size_t nSweeps){
  // Setup the random bits
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<double> CiNorm(basis.size(), 0.);
    for (size_t i=0; i<basis.size(); ++i){
      CiNorm[i] = std::norm(para.getCoeff(basis.getDetByIndex(i)));
    }
    std::discrete_distribution<> d(CiNorm.begin(), CiNorm.end());
    std::vector<int> count(basis.size(),0);  // holds count of each encountered number 
	  for(size_t i = 0; i<nSweeps; ++i){
      int r = d(gen);
      //count[d(gen)]++; 
      count[r]++; 
    }
    return count;
}

std::vector<int> genMarkovSample(Sampler &sampler, 
    Basis const& basis, size_t nSweeps){
	  coeffType cCoeff = coeffType();
	  detType cDet = detType();
	  double cWeight;
    //set nSweeps larger than the Hilbert space, so that each det could be sampled
    std::vector<int> count(basis.size(),0);  // holds count of each encountered number 
	  for(size_t i = 0; i<nSweeps; ++i){
		  sampler.iterate(cCoeff,cDet, cWeight, i);
      count[basis.getIndexByDet(cDet)]++;
	  }
	  sampler.updateBiases();
    return count;
}

void runMetropolisTest(Parametrization const& para,
    Basis const& basis, Sampler &sampler, size_t nSweeps){
  std::vector<int> directSample = genDirectSample(para, basis, nSweeps);
  std::vector<int> markovSample = genMarkovSample(sampler, basis, nSweeps);
  std::cout << "L2 Norm Closeness= " << testL2Norm(directSample, markovSample, nSweeps) << std::endl;
  std::cout << "L1 Norm Closeness= " << testL1Norm(directSample, markovSample, nSweeps) << std::endl;
}

int main(){
  // set up the system
  int numSites{6};
  int numStates{numSites*2};
  int numHidden(10);
  // then the Hamiltonian
  auto basis = generateDefaultBasis(numSites);
  auto modelHam = generateDefaultHubbard(numSites);
  size_t nSweeps = std::max(1000 * basis.size(), size_t(10000));
  // The cost function does not matter, we only need the para to get the coeffs
  DirectParametrization para(basis);
  // RBM para(numStates, numHidden);

  detType HF=basis.getDetByIndex(0);
  std::cout << "IntCast of HF det=" << verbatimCast(HF) << std::endl;
  // and the sampler
  // run 4 tests: one for the weighted Excitgen
  //WeightedExcitgen wEG(modelHam,HF);
  //MetropolisSampler<std::complex<double>> weightedSampler(wEG,HF, basis, para);

  for (int i(0); i<1000; ++i){
    // one for the RSHubbard
    //RSHubbardExcitgen hubbardEG{};
    //MetropolisSampler<std::complex<double>> hubSampler(hubbardEG,HF, basis, para);
    //std::cout << "Closeness for RSHubbardExcitgen: "<<  std::endl;
    //runMetropolisTest(para, basis, hubSampler, nSweeps);
    // one for RSHubbard via default initialization
    //MetropolisSampler<std::complex<double>> sampler(modelHam,HF, basis, para);
    //std::cout << "Closeness for RSHubbardExcitgen via default initialization: "<< std::endl;
    //runMetropolisTest(para, basis, sampler, nSweeps);

    ////// and one for the Uniform
    UniformExcitgen uniEG(HF);
    MetropolisSampler ugSampler(uniEG, HF, basis, para);


    std::cout << "Closeness for uniform EG: "<< std::endl;
    runMetropolisTest(para, basis, ugSampler, nSweeps);
  }
}


