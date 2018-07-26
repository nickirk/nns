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
  int numHidden(20);
  double trainRate(0.005);
  Basis basis(spinConfig);
  FermiHubbardHamiltonian modelHam(numStates);
  double U{4}, t{1};
  modelHam = generateFermiHubbard(numStates, U, t);
  vector<detType> list;
  //generate hamiltonian
  //AbInitioHamiltonian modelHam(numStates);
  //double U{2.}, t{-1};
  //string file_name = "FCIDUMP";
  //modelHam = readAbInitioHamiltonian(numStates, file_name);
  RBM rbm(numStates, numHidden);

  detType HF=basis.getDetByIndex(0);
  //EnergyEsMarkov eCF(modelHam);
  // EnergyEsMarkov cost funciton
  // works only with Markov Chain sampler.
  // Don't not use it for other samplers.
  //MetropolisSampler<VecCType> sampler(modelHam, basis, HF, rbm);
  //sampler.setNumDets(1000);
  RSHubbardExcitgen RSHG;
  //UniformExcitgen RSHG(HF);
  //WeightedExcitgen RSHG(modelHam,HF);
  MetropolisSampler<VecCType> sampler(RSHG, basis, HF,rbm);
  //ListGen<VecCType> sampler(RSHG, basis, HF,rbm,100);
  sampler.setNumDets(1000);
  //FullSampler<VecCType> sampler(modelHam, basis, HF, rbm);
  EnergyEs eCF(modelHam,-1);
  //Setup the trainer
  double energy{0.0};
  //AcceleratedGradientDescent<VecCType> sl(trainRate);
  ADAM<VecCType> sl(trainRate);
  //StochasticReconfiguration<VecCType> sl(rbm,trainRate);
  Trainer<VecCType> ev(rbm, sampler, sl, eCF,modelHam);
  ofstream myfile1;
  myfile1.open ("en1");
  for(int l(0); l<100000; ++l){
    trainRate = std::max(0.002*std::pow(0.99,l), 0.001);
    ev.train(trainRate);
    // get the new energy
    energy = ev.getE();
    auto states=ev.getState();
    for(size_t s=0; s<states.size(); ++s){
      cout << "C_" << s << "= " << states.coeff(s) << endl;
    }
    std::cout << "iteration=" << l << std::endl;
    std::cout << "trainRate=" << trainRate << std::endl;
    std::cout<<"Energy "<<energy<<std::endl;
    myfile1 << l << "  " << energy << std::endl;
    // update the list of determinants used in the sampler
    //sampler.diffuse(list,spinConfig);
  }
  myfile1.close();
}
