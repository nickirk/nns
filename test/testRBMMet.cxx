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
  int numHidden(20);
  double trainRate(0.0001);
  double U{4.}, t{-1};
  FermiHubbardHamiltonian modelHam(U,t,numStates);
  //AbInitioHamiltonian modelHam(numStates);
  Basis basis(spinConfig,modelHam);
  //double U{2.}, t{-1};
  //string file_name = "FCIDUMP_abinitio";
  //modelHam = readAbInitioHamiltonian(numStates, file_name);
  //vector<detType> list;
  RBM rbm(numStates, numHidden);

  detType HF=basis.getDetByIndex(0);
  RSHubbardExcitgen RSHG;
  //WeightedExcitgen weg(modelHam,HF);
  //ListGen<std::complex<double>> ugSampler(modelHam, basis, HF, rbm);
  //ListGen<std::complex<double>> ugSampler(RSHG, basis, HF, rbm);
  //ugSampler.setNumDets(100);
  EnergyEs eCF(modelHam, -1);
  MetropolisSampler ugSampler(RSHG, HF, basis, rbm);
  ugSampler.setNumDets(500);
  //sampler.diffuse(list,spinConfig);
  //Setup the trainer
  coeffType energy{0.0};
  //AcceleratedGradientDescent<std::complex<double>> sl(trainRate);
  ADAM sl(trainRate);
  //StochasticReconfiguration<std::complex<double>> sl(rbm,trainRate);
  Trainer ev(rbm, ugSampler, sl, eCF,modelHam);
  ofstream myfile1;
  myfile1.open ("en2");
  for(int l(0); l<50000; ++l){
    //trainRate *= exp(-0.0002);
    std::cout << "trainRate=" << trainRate << std::endl;
    ev.train(trainRate);
    // get the new energy
    energy = ev.getE();
    //auto states=ev.getState();
    //for(size_t s=0; s<states.size(); ++s){
    //  cout << "C_" << s << "= " << states.coeff(s) << endl;
    //}
    std::cout << "iteration=" << l << std::endl;
    std::cout<<"Energy "<<energy<<std::endl;
    myfile1 << l << "  " << energy << std::endl;
    // update the list of determinants used in the sampler
    //sampler.diffuse(list,spinConfig);
  }
  myfile1.close();
}
