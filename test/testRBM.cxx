#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "../src/NNWLib.hpp"

using namespace Eigen;
using namespace networkVMC;
using namespace std;

int main(){
  int numSites(6);
  int numStates(2*numSites);
  int spinUp(3);
  int spinDown(3);
  SpinConfig spinConfig(spinUp, spinDown, numStates);
  int numHidden(5);
  double trainRate(0.005);
  double U{2}, t{1};
  FermiHubbardHamiltonian modelHam(U,t,numSites);
  Basis basis(spinConfig,modelHam);
  vector<detType> list;
  //generate hamiltonian
  AbInitioHamiltonian modelHam(numStates);
  double U{2.}, t{-1};
  string file_name = "FCIDUMP";
  modelHam = readAbInitioHamiltonian(file_name);
  RBM rbm(numStates, numHidden);

  std::cout << "Initial vals of par=" << std::endl;
  std::cout << rbm.pars() << std::endl;

  detType HF=basis.getDetByIndex(0);
  //EnergyEsMarkov eCF(modelHam);
  // EnergyEsMarkov cost funciton
  // works only with Markov Chain sampler.
  // Don't not use it for other samplers.
  //MetropolisSampler<VecCType> sampler(modelHam, basis, HF, rbm);
  //sampler.setNumDets(1000);
  //RSHubbardExcitgen RSHG;
  //UniformExcitgen RSHG(HF);
  //WeightedExcitgen RSHG(modelHam,HF);
  MetropolisSampler<VecCType> sampler(RSHG, HF,rbm);
  //ListGen<VecCType> sampler(RSHG, basis, HF,rbm,100);
  //sampler.setNumDets(1000);
  FullSampler<VecCType> sampler(modelHam, basis, rbm);
  EnergyEs eCF(modelHam,-1);
  //Setup the trainer
  coeffType energy{0.0};
  //AcceleratedGradientDescent<VecCType> sl(trainRate);
  //ADAM<VecCType> sl(trainRate);
  StochasticReconfiguration<VecCType> sl(rbm,trainRate);
  Trainer<VecCType> ev(rbm, sampler, sl, eCF,modelHam);
  ofstream myfile1;
  myfile1.open ("en");
  for(int l(0); l<20000; ++l){
    trainRate = std::max(0.1*std::pow(0.995,l), 0.1);
    ev.train(trainRate);
    // get the new energy
    energy = ev.getE();
    auto states=ev.getState();
    for(size_t s=0; s<states.size(); ++s){
      cout << "C_" << s << "= " << states.coeff(s) << endl;
    }
    std::cout << "iteration=" << l << std::endl;
    std::cout << "trainRate=" << trainRate << std::endl;
    std::cout<<"Energy real="<< std::real(energy)<<std::endl;
    std::cout<<"Energy imag="<< std::imag(energy)<<std::endl;
    myfile1 << l << "  " << std::real(energy) << std::endl;
    // update the list of determinants used in the sampler
    //sampler.diffuse(list,spinConfig);
  }
  myfile1.close();
}
