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
  int numHidden(10);
  double trainRate(0.001);
  Basis basis(spinConfig);
  FermiHubbardHamiltonian modelHam(numStates);
  double U{4.}, t{-1};
  modelHam = generateFermiHubbard(numStates, U, t);
  vector<detType> list;
  RBM rbm(numStates, numHidden);

  detType HF=basis.getDetByIndex(0);
  //State target(std::vector<detType>(basis.getSize(),HF), 
  //    std::vector<coeffType>(basis.getSize(),coeffType(10.,1.)));
  EnergyEstimator eCF(modelHam);
  //NormCF eCF(target);
  //cout << "HF intCast=" << verbatimCast(HF) << endl;
  //list.push_back(HF);
  //ListGen<VecCType> sampler(modelHam, basis, HF, rbm, 200);
  FullSampler<VecCType> sampler(modelHam, basis, HF, rbm);
  //MetropolisSampler<VecCType> sampler(modelHam, basis, HF, rbm);
  //ListGen sampler(modelHam, basis, 100, HF,rbm);
  //sampler.diffuse(list,spinConfig);
  //Setup the trainer
  double energy{0.0};
  //AcceleratedGradientDescent<VecCType> sl(trainRate);
  ADAM<VecCType> sl(trainRate);
  Trainer<VecCType> ev(rbm, sampler, sl, eCF,modelHam);
  ofstream myfile1;
  myfile1.open ("en1");
  for(int l(0); l<1000; ++l){
    //trainRate *= exp(-0.0002);
    ev.train(trainRate);
    // get the new energy
    energy = ev.getE();
    auto states=ev.getState();
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
