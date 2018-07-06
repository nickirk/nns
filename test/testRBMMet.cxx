#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "../src/NNWLib.hpp"

using namespace Eigen;
using namespace networkVMC;
using namespace std;

int main(){
  int numSites(4);
  int numStates(2*numSites);
  int spinUp(2);
  int spinDown(2);
  SpinConfig spinConfig(spinUp, spinDown, numStates);
  int numHidden(10);
  double trainRate(0.0005);
  Basis basis(spinConfig);
  FermiHubbardHamiltonian modelHam(numStates);
  double U{4.}, t{-1};
  modelHam = generateFermiHubbard(numStates, U, t);
  vector<detType> list;
  RBM rbm(numStates, numHidden);

  detType HF=basis.getDetByIndex(0);
  RSHubbardExcitgen RSHG;
  MetropolisSampler<VecCType> ugSampler(RSHG, basis, HF, rbm);
  ugSampler.setNumDets(1000);
	EnergyEs eCF(modelHam, 10);
  //MetropolisSampler<VecCType> sampler(modelHam, basis, HF, rbm);
  //sampler.diffuse(list,spinConfig);
  //Setup the trainer
  double energy{0.0};
  //AcceleratedGradientDescent<VecCType> sl(trainRate);
  ADAM<VecCType> sl(trainRate);
  //StochasticGradientDescent<VecCType> sl(trainRate);
  Trainer<VecCType> ev(rbm, ugSampler, sl, eCF,modelHam);
  ofstream myfile1;
  myfile1.open ("en1");
  for(int l(0); l<5000; ++l){
    //trainRate *= exp(-0.0002);
    ev.train(trainRate);
    // get the new energy
    energy = ev.getE();
    //auto states=ev.getState();
    //for(size_t s=0; s<states.size(); ++s){
    //  cout << "C_" << s << "= " << states.coeff(s) << endl;
    //}
    std::cout << "iteration=" << l << std::endl;
    std::cout << "trainRate=" << trainRate << std::endl;
    std::cout<<"Energy "<<energy<<std::endl;
    myfile1 << l << "  " << energy << std::endl;
    // update the list of determinants used in the sampler
    //sampler.diffuse(list,spinConfig);
  }
  myfile1.close();
}
