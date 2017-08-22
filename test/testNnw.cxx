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
  int numStates(4);
  int numEle(2);
  int numHidden(30);
  bool readFromFile(false);
  Basis basis(numStates,numEle);
  Hamiltonian modelHam(0.5, 0.3, 0.2, basis);
  if (readFromFile){
    modelHam.readFromFile("modelHam.txt");  
  }
  else{
    modelHam.writeToFile("modelHam.txt");
  }
  cout << "Basis size= " << basis.getSize() << endl;
  cout << "Hamiltonian size= " << modelHam.getSize() << endl;
  cout << "print out Ham element" << endl;

  Eigen::MatrixXd H(Eigen::MatrixXd::Zero(basis.getSize(),basis.getSize()));
  int row, col;
  double value;
  for (int i(0); i<modelHam.getSparseSize(); ++i){
      modelHam.sparseAccess(i, row, col, value);
      H(row,col) = value; 
      cout <<"i= "<< i << endl ;
    }

  if(!readFromFile){
    ofstream myfile;
    myfile.open ("eigen.txt");
    cout << "writing eigen values to file" << endl;
    myfile << H.eigenvalues().real() << endl;

    //cout << H.eigenvalues().real() << endl;
    myfile.close();
  } 

  vector<int> size_NNW = {numStates, numHidden,1};
  vector<detType> list;
  //for (int i=0; i< basis.getSize(); ++i){
  //  list.push_back(basis.getDetByIndex(i));
  //  vector<int> pos=getOccupiedPositions(basis.getDetByIndex(i));
  //  for (int j=0; j<pos.size(); j++){
  //    cout << pos[j] << ",";
  //  }
  //  cout << endl;
  //}
  NeuralNetwork NNW(size_NNW, modelHam, basis);
  detType HF=basis.getDetByIndex(0);
  //list.push_back(HF); 
  int numDetsToTrain_ = 60;
  Sampler sampler(modelHam, basis, numDetsToTrain_, HF);
  ofstream myfile1;
  myfile1.open ("energy.txt");
  double aveEnergy(0.);
  double totalenergy(0.);
  double energy(0.);
  double lastEnergy(0.);
  double lastAveEnergy(0.);
  int sign(0);
  int lastSign(0);
  int count(0);
  int count1(0);
  //test sampler;
  sampler.generateList(list); 
  for (size_t i=0; i<list.size(); ++i){
      cout<<"intCast= " << verbatimCast(list[i])<<endl;
    }
  while (true){
    //list = NNW.train(list, 0.1); 
    //cout << "seeds size= " << list.size() << endl;
    sampler.generateList(list); 
    for (size_t i=0; i<list.size(); ++i){
      cout<<"intCast= " << verbatimCast(list[i])<<endl;
    }
    lastAveEnergy = aveEnergy;
    lastEnergy = energy;
    lastSign = sign;
    NNW.train(list, 0.08);
    energy = NNW.getEnergy();
    sign = NNW.getSign();
    count++;
    if (count1 < 50){
      totalenergy+=energy; 
      count1++;
      aveEnergy = totalenergy/count1;
    }
    else{
      totalenergy+=energy; 
      count1++;
      aveEnergy = totalenergy/count1;
      totalenergy=0.;
      count1=0;
    }
    cout << "Ave energy= " << aveEnergy<< endl;
    cout << "energy= " << energy<< endl;
    cout << "list size= " << list.size()<< endl;
    cout << "sign = " << sign<< endl;
    myfile1 << count << " " << energy << " " << aveEnergy << endl;
    

    cout << "percentage of allowed sign change= " << int(numDetsToTrain_*0.4) << endl;
    cout << "sign= " << sign << endl;
    cout << "lastSign= " << lastSign << endl;
    cout << "number of sign changes= " << abs(sign-lastSign) << endl;
    //cout << "abs(energy - lastEnergy)= " << abs(energy - lastEnergy) << endl;
    //cout << "fabs(energy - lastEnergy)= " << fabs(energy - lastEnergy) << endl;
    //cout << "abs(lastAveEnergy - aveEnergy)= " << abs(energy - lastEnergy) << endl;
    //while (abs(lastAveEnergy - aveEnergy) < 0.001 || abs(lastEnergy-energy) < 0.001){
    //while (abs(energy - lastEnergy) < 0.0001){
    //while (abs(sign-lastSign) > int(numDetsToTrain_*0.6) || abs(lastAveEnergy-aveEnergy) < 0.0005){
    while (false){
      list = NNW.train(list, 0.08);
      sampler.setReference(list);
      cout << "size of seeds= " << list.size() << endl;
      for (size_t i=0; i<list.size(); ++i){
        cout<<"Seeds intCast= " << verbatimCast(list[i])<<endl;
      }
      if (list.size() == numDetsToTrain_){ 
       numDetsToTrain_+=numDetsToTrain_;
       sampler.setNumStates(numDetsToTrain_);
      }
      energy = NNW.getEnergy();
      lastSign = sign;
      sign = NNW.getSign();
      sampler.generateList(list);
      cout << "sign = " << sign<< endl;
      //totalenergy+=energy; 
      count++;
      aveEnergy = energy/count;
    }
  }
  myfile1.close();
}
