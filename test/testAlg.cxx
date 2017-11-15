#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include "../src/Nnw.hpp"
#include "../src/Basis.hpp"
#include "../src/Determinant.hpp"
#include "../src/Hamiltonian.hpp"
#include "../src/Sampler.hpp"
#include "../src/EnergyEstimator.hpp"
#include "../src/EnergyCF.hpp"
using namespace Eigen;

using namespace std;
int main(){
  int numSites(8);
  int numStates(2*numSites);
  int numEle(8);
  int spinUp(8);
  int spinDown(0);
  vector<int> spinConfig{spinUp, spinDown, numStates};
  int numHidden(10*numSites);
  int numHidden1(3);
  vector<int> size_NNW = {numStates, numHidden,2};
  //cout << "input number of hidden neurons=";
  //cin >> numHidden;
  bool readFromFile{false};
  double trainRate(1.5);
  cout << "input training rate=";
  cin >> trainRate;
  Basis basis(spinConfig);
  Hamiltonian modelHam(numStates);
  double U{4.}, t{-1};
  modelHam = generateHubbard(numStates, U, t);
  
  cout << "Basis size= " << basis.getSize() << endl;
  //cout << "Hamiltonian size= " << modelHam.getSize() << endl;
  cout << "print out Ham element" << endl;
  vector<detType> list;
  ofstream detsIntcast; 

  detsIntcast.open("intCast.txt");
  for (int i=0; i< basis.getSize(); ++i){
    list.push_back(basis.getDetByIndex(i));
    detsIntcast << verbatimCast(basis.getDetByIndex(i)) << endl;
  }
  detsIntcast.close();
  std::cout<<"Listsize= "<<list.size()<<std::endl;
  //for(size_t i = 0; i< list.size(); ++i){
  //  std::cout<<"intCast= "<<verbatimCast(list[i])<<std::endl;
  //}
  EnergyEstimator eCF(modelHam);
  NeuralNetwork NNW(size_NNW, eCF);
  detType HF=basis.getDetByIndex(0);
  //cout << "HF intCast=" << verbatimCast(HF) << endl;
  //list.push_back(HF); 
  int numDetsToTrain_ = basis.getSize();
  cout << "numDetsToTrain= ";
  cin >> numDetsToTrain_;
  Sampler sampler(modelHam, basis, NNW, numDetsToTrain_, HF);
  ofstream myfile1;
  myfile1.open ("energy.txt");
  double aveEnergy(0.);
  double totalenergy(0.);
  double energy(0.);
  int sign(0);
  int lastSign(0);
  int count(0);
  int count1(0);
  std::vector<Eigen::VectorXd> coeffs;
  //test sampler;
  // set reference list Dets
  vector<detType> refDets;
  //refDets.push_back({1,0,0,1,1,0,0,1,1,0,0,1});
  //refDets.push_back({0,1,1,0,0,1,1,0,0,1,1,0});
  //sampler.setReference(refDets);
  //sampler.initialiseList(list, spinConfig); 
  //for (size_t i=0; i<list.size(); ++i){
  //  }
  double epsilon(0.3);
  double aveEnergyPrev(0.);
  double energyPrev(0.);
  int refSize(0);
  int listSize(0);
  vector<detType> listRef;
  vector<detType> listRefPrev;
  vector<detType> listRefTotal;
  while (true){
    NNW.train(list, trainRate, epsilon);
    sampler.diffuse(list, spinConfig); 
    listSize = list.size();
    energy = NNW.getEnergy();
    count++;
    double aveCount = 10;
    if(U<3.999){
      U+=0.001;
      modelHam = generateHubbard(numStates, U, t);
    }
    energyPrev = energy;
    if (epsilon > 0.05)
      epsilon -= 0.001;
    if(count%5 == 0){
      State states=NNW.getState();
      double normalizer=eCF.getNormalizer();
      ofstream outputC;
      outputC.open("coeff.txt");
      for(size_t s=0; s<states.size(); ++s){
        outputC << verbatimCast(states.getDet(s)) << " " << sqrt(norm(states.getCoeff(s)))/sqrt(normalizer) << endl; 
       //cout << "C_" << s << "= " << states.getCoeff(s) << endl;
      }
      cout << "normalizer=" << normalizer << endl;
      outputC.close();
      myfile1 << count << " " << energy << " " <<   " " << aveEnergy<< endl;
      int allowedNumChangeSign = int(basis.getSize()*0.1);
      cout << "energy= " << energy<< endl;
      cout << "list size= " << list.size()<< endl;
      cout << "U= " << U << endl;
      if (abs(sign-lastSign) > allowedNumChangeSign){
        trainRate*=0.90;
        cout << "trainRate=" << trainRate << endl;
      }
    }
  }
  myfile1.close();
}
