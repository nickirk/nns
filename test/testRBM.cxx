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
  int spinUp(3);
  int spinDown(3);
  int numHidden(24);
  double trainRate(0.005);
  double U{4}, t{1};
  int numStates = numSites*2; 
  FermiHubbardHamiltonian modelHam(U,t,numSites);
  //HeisenbergHamiltonian modelHam(-1.,false, numSites);
  //vector<detType> list;
  //generate hamiltonian
  //AbInitioHamiltonian modelHam(0);
  //string file_name = "FCIDUMP";
  //modelHam = readAbInitioHamiltonian(file_name,1);
  //numStates = modelHam.getNumOrbs();
  //std::cout << "numStates=" << numStates << std::endl;
  SpinConfig spinConfig(spinUp, spinDown,numStates);// numStates);
  Basis basis(spinConfig,modelHam);
  std::cout << "basis size=" << basis.size() << std::endl;
  //RBM rbm(numStates, numHidden);
  DirectParametrization<VecCType> rbm(basis);

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
  //LatticeExcitgen RSHG(modelHam);
  UniformExcitgen RSHG(HF);
  //WeightedExcitgen RSHG(modelHam,HF);
  MetropolisSampler<VecCType> sampler(RSHG, HF,basis, rbm);
  //ListGen<VecCType> sampler(RSHG, basis, HF,rbm,100);
  sampler.setNumDets(500);
  //FullSampler<VecCType> sampler(modelHam, basis, rbm);
  EnergyEs eCF(modelHam,-1);
  //Setup the trainer
  coeffType energy{0.0};
  //AcceleratedGradientDescent<VecCType> sl(trainRate);
  //ADAM<VecCType> sl(trainRate);
  StochasticReconfiguration<VecCType> sl(rbm,trainRate);
  Trainer<VecCType> ev(rbm, sampler, sl, eCF,modelHam);
  ofstream myfile1;
  myfile1.open ("en2");
  for(int l(0); l<200000; ++l){
    trainRate = std::max(0.5*std::pow(0.99,l), 0.1);
    ev.train(trainRate);
    // get the new energy
    energy = ev.getE();
    auto states=ev.getState();
    for(size_t s=0; s<states.size(); ++s){
      cout << "C_" << s << "= " << states.coeff(s)  << "  intCast="<< verbatimCast(states.det(s)) << endl;
    }
    std::cout << "iteration=" << l << std::endl;
    std::cout << "trainRate=" << trainRate << std::endl;
    std::cout<<"Energy real="<< std::real(energy)<<std::endl;
    std::cout<<"Energy imag="<< std::imag(energy)<<std::endl;
    myfile1 << l << "  " << std::real(energy) << std::endl;
    auto vals=rbm.pars();
    ofstream myfile2;
    myfile2.open ("val");
    for (size_t v(0); v<vals.size(); v++){
      myfile2 << v << "  " << vals(v).real() << "  " << vals(v).imag() << std::endl;
    }
    myfile2.close();
    // update the list of determinants used in the sampler
    //sampler.diffuse(list,spinConfig);
  }
  myfile1.close();
}
