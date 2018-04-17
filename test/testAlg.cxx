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
#include "../src/Trainer.hpp"
using namespace Eigen;

using namespace std;
int main(){
  int numSites(4);
  int numStates(2*numSites);
  int spinUp(2);
  int spinDown(2);
  
  vector<int> spinConfig{spinUp, spinDown, numStates};
  //int numHidden1(2*numSites);
  //cout << "input number of hidden neurons=";
  //cin >> numHidden;
  double trainRate(0.1);
  //cout << "input training rate=";
  //cin >> trainRate;
  //generate basis, the basis class constructor takes in the spin configurations.
  Basis basis(spinConfig);
  //generate hamiltonian
  FermionicHamiltonian modelHam(numStates);
  double U{4.}, t{-1};
  modelHam = generateFermiHubbard(numStates, U, t);
  
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
  EnergyEstimator eCF(modelHam);
  //Neural network takes in the size and the cost function.
  NeuralNetwork NNW(modelHam, eCF);
  NNW.constrInputLayer(numStates);
  //constrConvLayer(inputs, actFunc, lengthFilter, depthFilter, stride)
  cout << "Before constructing ConvLayer" << endl;
  NNW.constrConvLayer(NNW.getLayer(0)->getActs(), "Tanh", 6, 4, 2);
  NNW.constrDenseLayer(NNW.getLayer(1)->getActs(), "Tanh", 10);
  NNW.constrDenseLayer(NNW.getLayer(2)->getActs(), "Linear", 2);
  NNW.initialiseNetwork();
  cout << "After iniialise Network" << endl;
  // numDetsToTrain_ is the total number you want to keep in the list 
  // during the training process. By default it is the whole Hilbert space.
  // But for a stochastic diffuse process, much less is needed. One should 
  // test to see how many is suitable for different systems.
  int numDetsToTrain_ = basis.getSize();
  //cout << "numDetsToTrain= ";
  //cin >> numDetsToTrain_;
  int method(3);
  //cout << "Choose solver: 0, 1, 2, 3" << endl;
  //cout << "method 0: gradient descent" << endl;
  //cout << "method 1: stochastic reconfiguration (not working)" << endl;
  //cout << "method 2: RMSprop (not so stable)" << endl;
  //cout << "method 3: ADAM (default)" << endl;
  //cin >> method;
  string fileName("en");
  //cout << "energy file name=";
  //cin >> fileName;
  ListGen sampler(modelHam, basis, numDetsToTrain_, HF,NNW);
  //sampler.diffuse(list,spinConfig);
  //Setup the trainer
  Trainer ev(modelHam, NNW,sampler);
  ofstream myfile1;
  myfile1.open (fileName);
  double energy(0.);
  int count(0);
  std::vector<Eigen::VectorXd> coeffs;
  for(int l(0); l<10000; ++l){
    cout << "testAlg.cxx: iteration= " << l << endl;
    sampler.diffuse(list,spinConfig);
    ev.train(trainRate,method,l);

    // get the new energy
    energy = ev.getE();

    // update the list of determinants used in the sampler
    count++;
    if(count%1 == 0){
      std::vector<State > states=ev.getState();
      double normalizer=eCF.getNormalizer();
      ofstream outputC;
      outputC.open("coeff.txt");
      for(size_t s=0; s<states.size(); ++s){
        outputC << verbatimCast(states[s].det) << " " 
        << sqrt(norm(states[s].coeff))/sqrt(normalizer) << endl; 
        cout << "C_" << s << "= " << states[s].coeff << endl;
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
      cout << "U= " << U << endl;
    }
  }
  myfile1.close();
}
