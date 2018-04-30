#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "../src/CostFunctions/EnergyCF.hpp"
#include "../src/CostFunctions/EnergyEstimator.hpp"
#include "../src/Hamiltonian/FermionicHamiltonian.hpp"
#include "../src/Samplers/ListGen.hpp"
#include "../src/Network/Nnw.hpp"
#include "../src/HilbertSpace/Basis.hpp"
#include "../src/HilbertSpace/Determinant.hpp"
#include "../src/Samplers/MetropolisSampler.hpp"
#include "../src/Trainer.hpp"

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
  NeuralNetwork NNW(eCF);

  NNW.constrInputLayer(numStates);
  NNW.constrDenseLayer(NNW.getLayer(0)->getActs(), "Tanh", 10*numStates);
  NNW.constrDenseLayer(NNW.getLayer(1)->getActs(), "Linear", 2);
  NNW.initialiseNetwork();

  detType HF=basis.getDetByIndex(0);
  //cout << "HF intCast=" << verbatimCast(HF) << endl;
  //list.push_back(HF);
  MetropolisSampler sampler(modelHam, basis, HF,NNW);
  //ListGen sampler(modelHam, basis, 100, HF,NNW);
  //sampler.diffuse(list,spinConfig);
  //Setup the trainer
  double energy{0.0};
  Trainer ev(NNW,sampler);
  for(int l(0); l<100; ++l){
    ev.train(trainRate,3,l);
    // get the new energy
    energy = ev.getE();
    std::cout<<"Energy "<<energy<<std::endl;
    // update the list of determinants used in the sampler
    //sampler.diffuse(list,spinConfig);
  }
}
