#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "../src/NNWLib.hpp"

using namespace Eigen;
using namespace networkVMC;
using namespace std;

int main(){
  int numSites(20);
  int spinUp(10);
  int spinDown(10);
  int numHidden(10);
  double trainRate(0.005);
  double U{4}, t{1};
  int numStates = numSites; 
  //FermiHubbardHamiltonian modelHam(U,t,numSites);
  HeisenbergHamiltonian modelHam(-1.,false, numSites);
  //vector<detType> list;
  //generate hamiltonian
  //AbInitioHamiltonian modelHam(0);
  //double U{2.}, t{-1};
  //string file_name = "FCIDUMP";
  //modelHam = readAbInitioHamiltonian(file_name,1);
  //int numStates = modelHam.getNumOrbs();
  //std::cout << "numStates=" << numStates << std::endl;
  SpinConfig spinConfig(spinUp, spinDown,numStates);// numStates);
  Basis basis(spinConfig,modelHam);
  std::cout << "basis size=" << basis.size() << std::endl;
  RBM rbm(numStates, numHidden);
  //DirectParametrization<VecCType> rbm(basis);

  std::cout << "Initial vals of par=" << std::endl;
  std::cout << rbm.pars() << std::endl;
  std::cout << "num of par=" << std::endl;
  std::cout << rbm.getNumPars() << std::endl;

  detType HF=basis.getDetByIndex(0);
  //EnergyEsMarkov eCF(modelHam);
  // EnergyEsMarkov cost funciton
  // works only with Markov Chain sampler.
  // Don't not use it for other samplers.
  //MetropolisSampler<VecCType> sampler(modelHam, basis, HF, rbm);
  //sampler.setNumDets(1000);
  //RSHubbardExcitgen RSHG;
  LatticeExcitgen RSHG(modelHam);
  //UniformExcitgen RSHG(HF);
  //WeightedExcitgen RSHG(modelHam,HF);
  MetropolisSampler<VecCType> sampler(RSHG, HF,basis, rbm);
  //ListGen<VecCType> sampler(RSHG, basis, HF,rbm,100);
  sampler.setNumDets(1000);
  //FullSampler<VecCType> sampler(modelHam, basis, rbm);
  EnergyEs eCF(modelHam,-1);
  //Setup the trainer
  coeffType energy{0.0};
  //AcceleratedGradientDescent<VecCType> sl(trainRate);
  //ADAM<VecCType> sl(trainRate);
  StochasticReconfiguration<VecCType> sl(rbm,trainRate);
  Trainer<VecCType> ev(rbm, sampler, sl, eCF,modelHam);
  ofstream myfile1;
  myfile1.open ("en3");
  for(int l(0); l<20000; ++l){
    trainRate = std::max(0.01*std::pow(0.95,l), 0.001);
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
