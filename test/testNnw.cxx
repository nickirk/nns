#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include "../src/Nnw.hpp"
#include "../src/Basis.hpp"
#include "../src/Determinant.hpp"
#include "../src/FermionicHamiltonian.hpp"
#include "../src/ListGen.hpp"
#include "../src/EnergyEstimator.hpp"
#include "../src/EnergyCF.hpp"
#include "../src/MetropolisSampler.hpp"
#include "../src/Trainer.hpp"
using namespace Eigen;

using namespace std;
int main(){
  int numSites(6);
  int numStates(2*numSites);
  int spinUp(3);
  int spinDown(3);
  vector<int> spinConfig{spinUp, spinDown, numStates};
  int numHidden(40*numSites);
  vector<int> size_NNW = {numStates, numHidden, 2};
  double trainRate(0.1);
  Basis basis(spinConfig);
  FermionicHamiltonian modelHam(numStates);
  double U{4.}, t{-1};
  modelHam = generateFermiHubbard(numStates, U, t);
  vector<detType> list;
  EnergyEstimator eCF(modelHam);
  NeuralNetwork NNW(size_NNW, eCF);
  detType HF=basis.getDetByIndex(0);
  //cout << "HF intCast=" << verbatimCast(HF) << endl;
  //list.push_back(HF);
  int numDetsToTrain_{200};
  MetropolisSampler sampler(modelHam, basis, numDetsToTrain_, HF,NNW);
  //sampler.diffuse(list,spinConfig);
  //Setup the trainer
  double energy{0.0};
  Trainer ev(NNW,sampler);
  for(int l(0); l<10000; ++l){
    ev.train(trainRate,2);
    // get the new energy
    energy = ev.getE();
    std::cout<<"Energy "<<energy<<std::endl;
    // update the list of determinants used in the sampler
    //sampler.diffuse(list,spinConfig);
  }
}
