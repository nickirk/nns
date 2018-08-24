/*
 * testUmbrella.cxx
 *
 *  Created on: Aug 23, 2018
 *      Author: liao
 */
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
  //set up the system to study
  int numSites(4);
  int spinUp(2);
  int spinDown(2);
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

  // set up the frist RBM to capture the rough feature of 
  // the system
  //-----------------------------------------------------//
  // if read RBM from file
  bool readRBMFromFile= false;
  RBM<std::complex<double>, std::complex<double>> rbm(numStates, numHidden);
  if (readRBMFromFile){
    std::cout << "Reading RBM pars from file" << std::endl;
    rbm.readParsFromFile("Parameters.dat");
  }
  //DirectParametrization<std::complex<double>> rbm(basis);

  std::cout << "Initial vals of par=" << std::endl;
  std::cout << rbm.pars() << std::endl;
  std::cout << "num of par=" << std::endl;
  std::cout << rbm.getNumPars() << std::endl;

  // set up the excitation generator and MetropolisSampler
  // for the 1st round of training
  //-----------------------------------------------------//
  detType HF=basis.getDetByIndex(0);
  //UniformExcitgen RSHG(HF);
  LatticeExcitgen RSHG(modelHam);
  //WeightedExcitgen RSHG(modelHam,HF);
  MetropolisSampler<std::complex<double>,std::complex<double>> sampler(RSHG, HF,basis, rbm);
  //ListGen<std::complex<double>> sampler(RSHG, basis, HF,rbm,100);
  // set number of samplers in the Markov chain
  sampler.setNumDets(500);
  //FullSampler<std::complex<double>> sampler(modelHam, basis, rbm);

  // set up the cost function, solver and trainer
  //-----------------------------------------------------//
  EnergyEs<std::complex<double>,std::complex<double>> eCF(modelHam,-1);
  //Setup the trainer
  std::complex<double> energy{0.0};
  //AcceleratedGradientDescent<std::complex<double>> sl(trainRate);
  //ADAM<std::complex<double>, std::complex<double>> sl(trainRate);
  StochasticReconfiguration<std::complex<double>> sl(rbm,trainRate);
  Trainer<std::complex<double>, std::complex<double>> ev(rbm, sampler, sl, eCF,modelHam);
  ofstream myfile1;

  // open file to record the energy
  myfile1.open ("en");
  // decide if reading parameters from file, if not train a RBM as trial wf
  // Start training for a short time
  if (!readRBMFromFile){
    for(int l(0); l<300; ++l){
      trainRate = std::max(0.01*std::pow(0.999,l), 0.01);
      ev.train(trainRate);
      // get the new energy
      energy = ev.getE();
      std::cout << "iteration=" << l << std::endl;
      std::cout << "trainRate=" << trainRate << std::endl;
      std::cout<<"Energy real="<< std::real(energy)<<std::endl;
      std::cout<<"Energy imag="<< std::imag(energy)<<std::endl;
      myfile1 << l << "  " << std::real(energy) << std::endl;

      // write to file the sampled wavefunction and rbm wavefunction
      // this can slow down the simulation a lot. Use only in small calcs
      // Use it only to examine the behaviour of the algorithm
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
    // write Trained parameters to file after training for the next
    // round of training
    rbm.writeParsToFile("Parameters.dat");
    std::cout << "Trained Pars=\n" << rbm.pars() << std::endl;
  }

  // start the 2nd round of training, set up another
  // rbm and trialWf
  //-----------------------------------------------------//

  // set up another randomly initialised rbm
  RBM<std::complex<double>, std::complex<double>> rbm1(numStates, numHidden);


  // combine two para to a single one
	TrialWfPara<std::complex<double>> twfPara(rbm1,rbm);

  // also a new sampler and a new trainer
  UmbrellaSampler<> USampler(RSHG, HF,basis, twfPara);
  //StochasticReconfiguration<std::complex<double>> slU(twfPara,trainRate);
  ADAM<std::complex<double>, std::complex<double>> slU(trainRate);
  Trainer<std::complex<double>, std::complex<double>> evU(twfPara, USampler, slU, eCF,modelHam);
  // set up the training process for twfPara, similar to the first rbm
  
  for(int l(0); l<10000; ++l){
    trainRate = std::max(0.1*std::pow(0.999,l), 0.01);
    evU.train(trainRate);
    // get the new energy
    energy = evU.getE();
    std::cout << "iteration umbrella=" << l << std::endl;
    std::cout << "trainRate=" << trainRate << std::endl;
    std::cout<<"Energy real="<< std::real(energy)<<std::endl;
    std::cout<<"Energy imag="<< std::imag(energy)<<std::endl;
    myfile1 << l+300 << "  " << std::real(energy) << std::endl;

    // write to file the sampled wavefunction and rbm wavefunction
    // this can slow down the simulation a lot. Use only in small calcs
    // Use it only to examine the behaviour of the algorithm
    auto states=evU.getState();
    std::vector<int> intCasts(states.size());
    for(size_t s=0; s<states.size(); ++s){
      intCasts[s]=basis.getIndexByDet(states.det(s));
      //intCasts[s]=verbatimCast(states.det(s));
    }
    std::unordered_map<int, size_t> count;  // holds count of each encountered number 
    for (size_t i=0; i<intCasts.size(); i++)        
      count[intCasts[i]]++;      

    ofstream myfile3;
    myfile3.open ("SampleFunc");
    for (auto &e:count){
      myfile3 <<  e.first << " " << e.second/double(intCasts.size()) << std::endl;
    }
    myfile3.close();

    double normalisation(0.);
    for (size_t v(0); v<basis.size(); v++){
      normalisation+=std::norm(twfPara.getBaseCoeff(basis.getDetByIndex(v)));
    }
    ofstream myfile2;
    myfile2.open ("WaveFunc");
    for (size_t v(0); v<basis.size(); v++){
      detType det = basis.getDetByIndex(v);
      double amp = std::norm(twfPara.getBaseCoeff(det));
      std::cout << "norm(coeff(" << v << "))= " << amp << std::endl;
      myfile2 << v << "  " << amp/normalisation << std::endl;
    }
    myfile2.close();

    std::cout << "twfPara Base=\n" << twfPara.pars() << std::endl;
  }
  myfile1.close();

  //std::cout << "reading rbm pars from file..." << std::endl;
  //rbm.readParsFromFile("Parameters.dat");

  //rbm.readParsFromFile("Parameters.dat");
  //std::cout << "Read Pars=\n" << rbm.pars() << std::endl;
}
