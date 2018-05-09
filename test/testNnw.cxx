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
  int numHidden(40*numSites);
  vector<int> size_NNW = {numStates, numHidden, 2};
  double trainRate(0.1);
  Basis basis(spinConfig);
  FermionicHamiltonian modelHam(numStates);
  double U{4.}, t{-1};
  modelHam = generateFermiHubbard(numStates, U, t);
  vector<detType> list;
  EnergyEstimator eCF(modelHam);
  NeuralNetwork<VecType> NNW;

  NNW.constrInputLayer(numStates);
  NNW.constrDenseLayer(NNW.getLayer(0)->getActs(), "Tanh", 10*numStates);
  NNW.constrDenseLayer(NNW.getLayer(1)->getActs(), "Linear", 2);
  NNW.initialiseNetwork();

  detType HF=basis.getDetByIndex(0);
  //cout << "HF intCast=" << verbatimCast(HF) << endl;
  //list.push_back(HF);
  ListGen<> sampler(modelHam, basis, HF,NNW);
  //ListGen sampler(modelHam, basis, 100, HF,NNW);
  //sampler.diffuse(list,spinConfig);
  //Setup the trainer
  double energy{0.0};
  ADAM<VecType> sl(trainRate);
  Trainer<VecType> ev(NNW,sampler,sl,eCF);
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
