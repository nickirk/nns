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
  int numSites(8);
  int numStates(2*numSites);
  int spinUp(4);
  int spinDown(4);
  vector<int> spinConfig{spinUp, spinDown, numStates};
  int numHidden(40*numSites);
  int numHidden1(2*numSites);
  vector<int> size_NNW = {numStates, numHidden, 2};
  //cout << "input number of hidden neurons=";
  //cin >> numHidden;
  double trainRate(1.5);
  cout << "input training rate=";
  cin >> trainRate;
  Basis basis(spinConfig);
  FermionicHamiltonian modelHam(numStates);
  double U{4.}, t{-1};
  modelHam = generateFermiHubbard(numStates, U, t);
  
  cout << "Basis size= " << basis.getSize() << endl;
  //cout << "Hamiltonian size= " << modelHam.getSize() << endl;
  vector<detType> list;
  ofstream detsIntcast; 

  detsIntcast.open("intCast.txt");
  for (int i=0; i< basis.getSize(); ++i){
    list.push_back(basis.getDetByIndex(i));
    detsIntcast << verbatimCast(basis.getDetByIndex(i)) << endl;
  }
  detsIntcast.close();
  std::cout<<"Listsize= "<<list.size()<<std::endl;
  EnergyEstimator eCF(modelHam);
  NeuralNetwork NNW(size_NNW, eCF);
  detType HF=basis.getDetByIndex(0);
  //cout << "HF intCast=" << verbatimCast(HF) << endl;
  //list.push_back(HF); 
  int numDetsToTrain_ = basis.getSize();
  cout << "numDetsToTrain= ";
  cin >> numDetsToTrain_;
  int method(3);
  cout << "Choose solver: 0, 1, 2, 3" << endl;
  cout << "method 0: gradient descent" << endl;
  cout << "method 1: stochastic reconfiguration (not working)" << endl;
  cout << "method 2: RMSprop (not so stable)" << endl;
  cout << "method 3: ADAM (default)" << endl;
  cin >> method;
  string fileName;
  cout << "energy file name=";
  cin >> fileName;
  ListGen sampler(modelHam, basis, numDetsToTrain_, HF,NNW);
  //sampler.diffuse(list,spinConfig);
  //Setup the trainer
  Trainer ev(NNW,sampler);
  ofstream myfile1;
  myfile1.open (fileName);
  double energy(0.);
  int lastSign(0);
  int count(0);
  int count1(0);
  std::vector<Eigen::VectorXd> coeffs;
  double epsilon(0.3);
  int listSize(0);
  vector<detType> listRef;
  vector<detType> listRefPrev;
  vector<detType> listRefTotal;
  for(int l(0); l<10000; ++l){
    ev.train(trainRate,method,l);
    listSize = list.size();

    // get the new energy
    energy = ev.getE();

    // update the list of determinants used in the sampler
    sampler.diffuse(list,spinConfig);
    count++;
    if(U<3.999){
      U+=0.001;
      modelHam = generateFermiHubbard(numStates, U, t);
    }
    if (epsilon > 0.05)
      epsilon -= 0.001;
    if(count%1 == 0){
      std::vector<State > states=ev.getState();
      double normalizer=eCF.getNormalizer();
      ofstream outputC;
      outputC.open("coeff.txt");
      for(size_t s=0; s<states.size(); ++s){
        outputC << verbatimCast(states.getDet(s)) << " " << sqrt(norm(states.getCoeff(s)))/sqrt(normalizer) << endl; 
        cout << "C_" << s << "= " << states.getCoeff(s) << endl;
      }
      outputC.close();
      cout << "normalizer=" << normalizer << endl;
      myfile1 << count << " " << energy << endl;
      ofstream bias0;
      bias0.open ("bias0.txt");
      bias0 << NNW.getBiases(0) << endl;
      ofstream bias1;
      bias1.open ("bias1.txt");
      bias1 << NNW.getBiases(1) << endl;
      int allowedNumChangeSign = int(basis.getSize()*0.1);
      cout << "energy= " << energy<< endl;
      cout << "list size= " << list.size()<< endl;
      //cout << "weight 0=" << endl;
      //cout << NNW.getWeights(0) << endl;
      cout << "U= " << U << endl;
    }
  }
  myfile1.close();
}
