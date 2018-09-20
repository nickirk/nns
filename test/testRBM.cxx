#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <unordered_map> //for unordered_map
#include "../src/CostFunctions/EnergyEs.hpp"
#include "../src/NNWLib.hpp"

using namespace Eigen;
using namespace networkVMC;
using namespace std;

int main(){
  int numSites(6);
  int spinUp(3);
  int spinDown(3);
  int numHidden(10);
  double trainRate(0.005);
  double U{4}, t{-1};
  int numStates = numSites*2;

  FermiHubbardHamiltonian modelHam(U,t,numSites);
  //HeisenbergHamiltonian modelHam(-1.,true, numSites);
  //vector<detType> list;
  //generate hamiltonian
  //AbInitioHamiltonian modelHam(0);
  //std::string file_name = "FCIDUMP";
  //modelHam = readAbInitioHamiltonian(file_name,1);
  //numStates = modelHam.getNumOrbs();
  //std::cout << "numStates=" << numStates << std::endl;
  SpinConfig spinConfig(spinUp, spinDown,numStates);// numStates);
  Basis basis(spinConfig,modelHam);
  std::cout << "basis size=" << basis.size() << std::endl;
  //RBM<std::complex<double>, std::complex<double>> rbm(numStates, numHidden);
  DirectParametrization rbm(basis);

  std::cout << "Initial vals of par=" << std::endl;
  std::cout << rbm.pars() << std::endl;
  std::cout << "num of par=" << std::endl;
  std::cout << rbm.getNumPars() << std::endl;

  detType HF=basis.getDetByIndex(0);
  UniformExcitgen RSHG(HF);

  //LatticeExcitgen RSHG(modelHam);
  //WeightedExcitgen RSHG(modelHam,HF);
  MetropolisSampler sampler(RSHG, HF,basis, rbm);
  //ListGen sampler(RSHG, basis, HF,rbm,100);
  sampler.setNumDets(500);
  //FullSampler sampler(modelHam, basis, rbm);

  EnergyEs eCF(modelHam,-1);
  //Setup the trainer
  std::complex<double> energy{0.0};
  //AcceleratedGradientDescent sl(trainRate);
  //ADAM sl(trainRate);
  StochasticReconfiguration sl(rbm,trainRate);
  Trainer ev(rbm, sampler, sl, eCF,modelHam);
  ofstream myfile1;
  myfile1.open ("en");
  //std::cout << "reading rbm pars from file..." << std::endl;
  //rbm.readParsFromFile("Parameters.dat");
  for(int l(0); l<10000; ++l){
    trainRate = std::max(0.01*std::pow(0.999,l), 0.01);
    ev.train(trainRate);
    // get the new energy
    energy = ev.getE();
    std::cout << "iteration=" << l << std::endl;
    std::cout << "trainRate=" << trainRate << std::endl;
    std::cout<<"Energy real="<< std::real(energy)<<std::endl;
    std::cout<<"Energy imag="<< std::imag(energy)<<std::endl;
    myfile1 << l << "  " << std::real(energy) << std::endl;
    // update the list of determinants used in the sampler
    //sampler.diffuse(list,spinConfig);
    //Every 10 iterations write Pars to file
    auto states=ev.getState();
    std::vector<int> intCasts(states.size());
    for(size_t s=0; s<states.size(); ++s){
      intCasts[s]=basis.getIndexByDet(states.det(s));
      //intCasts[s]=verbatimCast(states.det(s));
    }
    std::unordered_map<int, size_t> count;  // holds count of each encountered number 
    for (size_t i=0; i<intCasts.size(); i++)        
      count[intCasts[i]]++;      

    double normalisation(0.);
    for (size_t v(0); v<basis.size(); v++){
      normalisation+=std::norm(rbm.getCoeff(basis.getDetByIndex(v)));
    }
    ofstream myfile3;
    myfile3.open ("SampleFunc");
    for (auto &e:count){
      myfile3 <<  e.first << " " << e.second/double(intCasts.size()) << std::endl;
    }
    myfile3.close();
    ofstream myfile2;
    myfile2.open ("WaveFunc");
    for (size_t v(0); v<basis.size(); v++){
      detType det = basis.getDetByIndex(v);
      double amp = std::norm(rbm.getCoeff(det));
      myfile2 << v << "  " << amp/normalisation << std::endl;
    }
    myfile2.close();
  }
  rbm.writeParsToFile("Parameters.dat");
  std::cout << "Trained Pars=\n" << rbm.pars() << std::endl;
  //rbm.readParsFromFile("Parameters.dat");
  //std::cout << "Read Pars=\n" << rbm.pars() << std::endl;
  myfile1.close();
}
