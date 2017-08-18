#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include "../src/Nnw.hpp"
#include "../src/Basis.hpp"
#include "../src/Determinant.hpp"
#include "../src/Hamiltonian.hpp"
#include "../src/Sampler.hpp"
using namespace std;

using namespace std;
int main(){
  int numStates(10);
  int numEle(3);
  Basis basis(numStates,numEle);
  Hamiltonian modelHam(0.5, 0.3, 0.2, basis);
  cout << "Basis size= " << basis.getSize() << endl;
  cout << "Hamiltonian size= " << modelHam.getSize() << endl;
  cout << "print out Ham element" << endl;
  Eigen::MatrixXd H(Eigen::MatrixXd::Zero(basis.getSize(),basis.getSize()));
  for (int i(0); i<modelHam.getSize(); ++i){
    for (int j(0); j< modelHam.getSize(); ++j){
      H(i,j) = modelHam(i,j); 
      cout << modelHam(i,j) << ", " ;
    }
    cout << endl;
  }
  ofstream myfile;
  myfile.open ("eigen.txt");
  cout << H.eigenvalues().real() << endl;
  myfile.close();
  vector<int> size_NNW = {numStates, 30,1};
  vector<detType> list;
  //for (int i=0; i<1; ++i){
  //  list.push_back(basis.getDetByIndex(i));
  //  vector<int> pos=getOccupiedPositions(basis.getDetByIndex(i));
  //  for (int j=0; j<pos.size(); j++){
  //    cout << pos[j] << ",";
  //  }
  //  cout << endl;
  //}
  NeuralNetwork NNW(size_NNW, modelHam, basis);
  detType HF=basis.getDetByIndex(0);
 
  int numDetsToTrain_ = 30;
  Sampler sampler(modelHam, basis, numDetsToTrain_, HF);
  sampler.generateList(list); 
  ofstream myfile1;
  myfile1.open ("energy.txt");
  double aveEnergy(0.);
  double totalenergy(0.);
  double energy(0.);
  double lastEnergy(0.);
  double lastAveEnergy(0.);
  int count(0);
  int count1(0);
  while (true){
    //list = NNW.train(list, 0.1); 
    //cout << "seeds size= " << list.size() << endl;
    //sampler.generateList(list); 
    lastAveEnergy = aveEnergy;
    lastEnergy = energy;
    NNW.train(list, 0.01);
    energy = NNW.getEnergy();
    count++;
    if (count1 < 300){
      totalenergy+=energy; 
      count1++;
      aveEnergy = totalenergy/count1;
    }
    else{
      totalenergy=0.;
      count1=0;
    }
    cout << "Ave energy= " << aveEnergy<< endl;
    cout << "energy= " << energy<< endl;
    cout << "list size= " << list.size()<< endl;
    myfile1 << count << " " << energy << " " << aveEnergy << endl;
    

    //while (abs(lastAveEnergy - aveEnergy) < 0.001 || abs(lastEnergy-energy) < 0.001){
    while (abs(energy - lastEnergy) < 0.0001){
      list = NNW.train(list, 0.08);
      cout << "size of seeds= " << list.size() << endl;
      if (list.size() == numDetsToTrain_){ 
       numDetsToTrain_=int(numDetsToTrain_+=1);
       sampler.setNumStates(numDetsToTrain_);
      }
      energy = NNW.getEnergy();
      sampler.generateList(list);
      //totalenergy+=energy; 
      count++;
      aveEnergy = energy/count;
    }
  }
  myfile1.close();
}
