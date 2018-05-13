#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "../src/NNWLib.hpp"
#include "../src/Trainer.cxx"

using namespace Eigen;
using namespace networkVMC;
using namespace std;

int main(){
  int numSites(6);
  int numStates(2*numSites);
  int spinUp(3);
  int spinDown(3);
  SpinConfig spinConfig(spinUp, spinDown, numStates);
  int numHidden(2*numSites);
  double trainRate(0.01);
  Basis basis(spinConfig);
  FermionicHamiltonian modelHam(numStates);
  double U{4.}, t{-1};
  modelHam = generateFermiHubbard(numStates, U, t);
  vector<detType> list;
  EnergyEstimator eCF(modelHam);
  RBM rbm(numStates, numHidden);

  detType HF=basis.getDetByIndex(0);
  //cout << "HF intCast=" << verbatimCast(HF) << endl;
  //list.push_back(HF);
  ListGen<VecCType> sampler(modelHam, basis, HF, rbm);
  //ListGen sampler(modelHam, basis, 100, HF,rbm);
  //sampler.diffuse(list,spinConfig);
  //Setup the trainer
  double energy{0.0};
  ADAM<VecCType> sl(trainRate);
  Trainer<VecCType> ev(rbm, sampler, sl, eCF);
  for(int l(0); l<100; ++l){
    ev.train();
    // get the new energy
    energy = ev.getE();
    std::cout << "iteration=" << l << std::endl;
    std::cout<<"Energy "<<energy<<std::endl;
    // update the list of determinants used in the sampler
    //sampler.diffuse(list,spinConfig);
  }
}
