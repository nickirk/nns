#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <Eigen/Dense>

#include "../src/NNWLib.hpp"
#include "defaultSystem.hpp"
using namespace Eigen;
using namespace networkVMC;

using namespace std;
int main(){
  int numSites = 6;
  int numStates = 2*numSites;
  int spinUp = 3;
  int spinDown = 3;
  int U = 4;
  int t = 1;
  //int numHidden1(2*numSites);
  //cout << "input number of hidden neurons=";
  //cin >> numHidden;
  double trainRate = 0.01;
  //cout << "input training rate=";
  //cin >> trainRate;
  //generate basis, the basis class constructor takes in the spin configurations.
  SpinConfig spinConfig{spinUp, spinDown, numStates};
  int numHidden(10);
  // which random excitation generator should be tested
  //generate basis, the basis class constructor takes in the spin configurations.
  FermionBasis basis(spinConfig);
  //generate hamiltonian
  AbInitioHamiltonian modelHam(numStates);
  ////double U{2.}, t{-1};
  string file_name = "FCIDUMP";
  modelHam = readAbInitioHamiltonian(file_name, numStates);
  //generate hamiltonian
  //FermiHubbardHamiltonian modelHam(numStates);
  //auto modelHam = generateFermiHubbard(numStates, U, t);
  cout << "Basis size= " << basis.getSize() << endl;
  //cout << "Hamiltonian size= " << modelHam.getSize() << endl;
  vector<detType> list;
  ofstream detsIntcast; 

  detsIntcast.open("intCast.txt");
  //in the diffuse scheme, we start with all the determinants in the basis
  //but this is not necessary. Because in the end we will diffuse them and
  //add random determinants.
  //IntCast just to cast the determinant into a integer number so that 
  //we can see if anything is going wrong with the basis.
  auto HF=basis.getDetByIndex(0);
  for (int i=0; i< basis.getSize(); ++i){
    list.push_back(basis.getDetByIndex(i));
    detsIntcast << verbatimCast(basis.getDetByIndex(i)) << endl;
  }
  detsIntcast.close();
  std::cout<<"Listsize= "<<list.size()<<std::endl;
  //Neural network takes in the size and the cost function.
  NeuralNetwork<> NNW;
  NNW.constrInputLayer(numStates);
  cout << "Before constructing ConvLayer" << endl;
  NNW.constrDenseLayer(NNW.getLayer(0)->getActs(), "Tanh", numHidden);
  NNW.constrDenseLayer(NNW.getLayer(1)->getActs(), "Linear", 2);
  NNW.initialiseNetwork();
  EnergyEs eCF(modelHam,-1);
  //WeightedExcitgen RSHG(modelHam, HF);
  //RSHubbardExcitgen RSHG;
  //UniformExcitgen RSHG(HF);
  //MetropolisSampler<VecType> sampler(RSHG, basis, HF,NNW);
  //sampler.setNumDets(1000);
  //ListGen<VecType> sampler(RSHG, basis, HF,NNW,1000);
  FullSampler<> sampler(modelHam, basis,NNW);
  ADAM<VecType> sl(trainRate);
  //AcceleratedGradientDescent<> sl(trainRate);
  //sampler.diffuse(list,spinConfig);
  //Setup the trainer
  Trainer<VecType> ev(NNW,sampler,sl,eCF,modelHam);
  string fileName("en");
  ofstream myfile1;
  myfile1.open (fileName);
  double energy(0.);
  int count(0);
  std::vector<Eigen::VectorXd> coeffs;
  for(int l(0); l<200000; ++l){
    cout << "testAlg.cxx: iteration= " << l << endl;
    //sampler.diffuse(list);
    //trainRate*=0.999;
    cout << "trainRate=" << trainRate << endl;
    ev.train(trainRate);

    // get the new energy
    energy = ev.getE();

    // update the list of determinants used in the sampler
    count++;
    if(count%1 == 0){
      auto states=ev.getState();
      // A horrible construct to get the normalizer of the trainer's cost function copy
      double normalizer=eCF.getNormalizer();
      ofstream outputC;
      outputC.open("coeff.txt");
      for(size_t s=0; s<states.size(); ++s){
        outputC << verbatimCast(states.det(s)) << " "
        << sqrt(norm(states.coeff(s)))/sqrt(normalizer) << endl;
        cout << "C_" << s << "= " << states.coeff(s) << endl;
      }
      outputC.close();
      cout << "normalizer=" << normalizer << endl;
      myfile1 << count << " " << energy << endl;
      /*
      ofstream bias0;
      bias0.open ("bias0.txt");
      bias0 << NNW.getBiases(0) << endl;
      ofstream bias1;
      bias1.open ("bias1.txt");
      bias1 << NNW.getBiases(1) << endl;
      */

      cout << "energy= " << energy<< endl;
      cout << "list size= " << list.size()<< endl;
      //cout << "weight 0=" << endl;
      //cout << NNW.getWeights(0) << endl;
    }
  }
  myfile1.close();
}
