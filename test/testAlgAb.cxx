#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include "../src/Nnw.hpp"
#include "../src/Basis.hpp"
#include "../src/Determinant.hpp"
#include "../src/FermionicHamiltonian.hpp"
#include "../src/AbInitioHamiltonian.hpp"
#include "../src/Sampler.hpp"
#include "../src/EnergyEstimator.hpp"
#include "../src/EnergyCF.hpp"
using namespace Eigen;

using namespace std;
int main(){
  int numSites(6);
  int numStates(2*numSites);
  int spinUp(3);
  int spinDown(3);
  
  vector<int> spinConfig{spinUp, spinDown, numStates};
  int numHidden(10*numSites);
  //int numHidden1(2*numSites);
  vector<int> size_NNW = {numStates, numHidden, 2};
  //cout << "input number of hidden neurons=";
  //cin >> numHidden;
  bool readFromFile{false};
  double trainRate(1.5);
  cout << "input training rate=";
  cin >> trainRate;
  //generate basis, the basis class constructor takes in the spin configurations.
  Basis basis(spinConfig);
  //generate hamiltonian
  AbInitioHamiltonian modelHam(numStates);
  string fileName="../run/FCIDUMP";
  modelHam = readAbInitioHamiltonian(numStates, fileName);
  cout << "Basis size= " << basis.getSize() << endl;
  //cout << "Hamiltonian size= " << modelHam.getSize() << endl;
  cout << "print out Ham element" << endl;
  vector<detType> list;
  ofstream detsIntcast; 

  detsIntcast.open("intCast.txt");
  //in the diffuse scheme, we start with all the determinants in the basis
  //but this is not necessary. Because in the end we will diffuse them and
  //add random determinants.
  //IntCast just to cast the determinant into a integer number so that 
  //we can see if anything is going wrong with the basis.
  for (int i=0; i< basis.getSize(); ++i){
    list.push_back(basis.getDetByIndex(i));
    detsIntcast << verbatimCast(basis.getDetByIndex(i)) << endl;
  }
  detsIntcast.close();
  std::cout<<"Listsize= "<<list.size()<<std::endl;
  //for(size_t i = 0; i< list.size(); ++i){
  //  std::cout<<"intCast= "<<verbatimCast(list[i])<<std::endl;
  //}
  //---------------------------------------------------//

  //define a energy estimator. We can explore more cost functions than just the 
  //energy in the future. That's why we make it a child class of the cost function.
  EnergyEstimator eCF(modelHam);
  //Neural network takes in the size and the cost function.
  NeuralNetwork NNW(modelHam, size_NNW, eCF);
  // numDetsToTrain_ is the total number you want to keep in the list 
  // during the training process. By default it is the whole Hilbert space.
  // But for a stochastic diffuse process, much less is needed. One should 
  // test to see how many is suitable for different systems.
  int numDetsToTrain_ = basis.getSize();
  cout << "numDetsToTrain= ";
  cin >> numDetsToTrain_;
  string energyFile;
  cout << "energy file name=";
  cin >> energyFile;

  //currently, the diffuse process is implemented within the sampler class.
  detType HF=list[0];
  Sampler sampler(modelHam, basis, NNW, numDetsToTrain_, HF);
  //write energy into this file so that we can plot the energy with the plotEn.sh script
  ofstream myfile1;
  myfile1.open (energyFile);
  double energy(0.);
  int count(0);
  std::vector<Eigen::VectorXd> coeffs;
  //test sampler;
  // set reference list Dets
  // there might be unused variables.
  vector<detType> refDets;
  double epsilon(0.3);
  double energyPrev(0.);
  int listSize(0);
  vector<detType> listRef;
  vector<detType> listRefPrev;
  vector<detType> listRefTotal;
  for(int l(0); l<10000; ++l){
    NNW.train(list, trainRate, l);
    sampler.diffuse(list, spinConfig); 
    listSize = list.size();
    energy = NNW.getEnergy();
    count++;
    energyPrev = energy;
    if (epsilon > 0.05)
      epsilon -= 0.001;
    if(count%1 == 0){
      State states=NNW.getState();
      double normalizer=eCF.getNormalizer();
      ofstream outputC;
      outputC.open("coeff.txt");
      for(size_t s=0; s<states.size(); ++s){
        outputC << verbatimCast(states.getDet(s)) << " " << sqrt(norm(states.getCoeff(s)))/sqrt(normalizer) << endl; 
       // cout << "C_" << s << "= " << states.getCoeff(s) << endl;
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
    }
  }
  myfile1.close();
}
