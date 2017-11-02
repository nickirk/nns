#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include "../src/Nnw.hpp"
#include "../src/Basis.hpp"
#include "../src/Determinant.hpp"
#include "../src/Hamiltonian.hpp"
#include "../src/Sampler.hpp"
#include "../src/EnergyCF.hpp"
using namespace Eigen;

using namespace std;
int main(){
  int numSites(6);
  int numStates(2*numSites);
  int numEle(6);
  int numHidden(8*numSites);
  int numHidden1(3);
  vector<int> size_NNW = {numStates, numHidden, 10, 2};
  //cout << "input number of hidden neurons=";
  //cin >> numHidden;
  bool readFromFile{false};
  double trainRate(0.005);
  cout << "input training rate=";
  cin >> trainRate;
  Basis basis(numStates,numEle);
  Hamiltonian modelHam(numStates);
  double U{4}, t{-1};
  modelHam = generateHubbard(numStates, U, t);
  
  cout << "Basis size= " << basis.getSize() << endl;
  //cout << "Hamiltonian size= " << modelHam.getSize() << endl;
  cout << "print out Ham element" << endl;
  /*
  Eigen::MatrixXd H(Eigen::MatrixXd::Zero(basis.getSize(),basis.getSize()));
  double sumRow(0.);
  for (int i(0); i<basis.getSize(); ++i){
          sumRow =0.;
          for(int j(0); j<basis.getSize(); ++j){
          H(i,j) = modelHam(basis.getDetByIndex(i),basis.getDetByIndex(j));
          sumRow += H(i,j); 
          }
  }
  Eigen::EigenSolver<MatrixXd> eSolver(H);
  VectorXd eVals = eSolver.eigenvalues().real();
  double eMin = 0.0;
  int pos = 0;
  for(int i = 0; i < basis.getSize(); ++i){
          if((eVals(i)-eMin)<1e-8){
              cout << "eigenVal="  << eVals(i) << " " << i << endl;
        	  pos = i;
        	  eMin = eVals(i);
          }
  }
  VectorXcd eVector = eSolver.eigenvectors().col(pos);
  //cout << eVector << endl;

  if(!readFromFile){
    ofstream myfile("eigen.txt");
    //myfile.open ("eigen.txt");
    cout << "writing eigen values to file" << endl;
    myfile << eVals <<  endl;

    myfile.close();
  } 
  if(!readFromFile){
    ofstream myfilevec;
    myfilevec.open ("eigenvec.txt");
    cout << "writing eigen vec to file" << endl;
    for (int i=0; i<basis.getSize(); ++i){
      myfilevec << verbatimCast(basis.getDetByIndex(i)) << " " << eVector(i).real() << endl;
    }
  }
 */
  vector<detType> list;
  ofstream detsIntcast; 
  detsIntcast.open("intCast.txt");
  for (int i=0; i< basis.getSize(); ++i){
    list.push_back(basis.getDetByIndex(i));
    detsIntcast << verbatimCast(basis.getDetByIndex(i)) << endl;
  }
  detsIntcast.close();
  std::cout<<"Listsize= "<<list.size()<<std::endl;
  for(size_t i = 0; i< list.size(); ++i){
    std::cout<<"intCast= "<<verbatimCast(list[i])<<std::endl;
  }
  EnergyCF eCF(modelHam);
  NeuralNetwork NNW(size_NNW, eCF);
  detType HF=basis.getDetByIndex(0);
  //list.push_back(HF); 
  int numDetsToTrain_ = basis.getSize();
  cout << "numDetsToTrain= ";
  cin >> numDetsToTrain_;
  Sampler sampler(modelHam, basis, numDetsToTrain_, HF);
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
  refDets.push_back({1,0,0,1,1,0,0,1,1,0,0,1});
  refDets.push_back({0,1,1,0,0,1,1,0,0,1,1,0});
  sampler.setReference(refDets);
  sampler.generateList(list); 
  sampler.removeDuplicate(list);
  //for (size_t i=0; i<list.size(); ++i){
  //    cout<<"intCast= " << verbatimCast(list[i])<<endl;
  //  }
  double maxEnergy(0.);
  double energySquare(0.);
  double variance(0.);
  double variancePrev(0.);
  double sampleEnergy(0.);
  double epsilon(0.05);
  double aveEnergyPrev(0.);
  double energyPrev(0.);
  int refSize(0);
  int listSize(0);
  while (true){
    //list = NNW.train(list, 0.1); 
    //cout << "seeds size= " << list.size() << endl;
    //for (size_t i=0; i<list.size(); ++i){
    //  cout<<"intCast= " << verbatimCast(list[i])<<endl;
    //}
    lastSign = sign;
    NNW.train(list, trainRate, epsilon);
    //list=NNW.train(list, trainRate, epsilon);
    //sampler.setReference(list);
    refSize = list.size();
    cout << "Ref list size= " << list.size()<< endl;
    //for (size_t i=0; i<list.size(); ++i){
    //  cout<<"Ref intCast= " << verbatimCast(list[i])<<endl;
    //}
    sampler.generateList(list); 
    sampler.removeDuplicate(list);
    cout << "New list size= " << list.size()<< endl;
    listSize = list.size();
    energy = NNW.getEnergy();
    count++;
    double aveCount = 1000;
    if (count1 < aveCount){
      totalenergy+=energy; 
      energySquare += pow(energy,2);
      count1++;
      if (count1 == aveCount-1){
        //if (maxEnergy - aveEnergy > 1) trainRate*=0.5;
       }
      //trainRate+=0.001;
    }
    else{
      list=NNW.train(list, trainRate, epsilon);
      sampler.setReference(list);
      sampler.generateList(list); 
      sampler.removeDuplicate(list);
      variancePrev = variance;
      aveEnergy = totalenergy/double(count1);
      if (fabs(aveEnergy - aveEnergyPrev) < 0.05) {
        //numDetsToTrain_ += 5;
        //sampler.setNumStates(numDetsToTrain_+1);
      }
      //if((refSize*1.0)/listSize < 0.5) {
      //  epsilon *=0.8;
      //  trainRate *= 0.95;
      //}
      //else trainRate *=1.05;
      totalenergy=0.;
      //count1 = 0.;
      energyPrev = energy;
      cout << "epsilon= " << epsilon << endl;
    }
    if(count%1 == 0){
    //cout << "sign = " << sign<< endl;
    State states=NNW.getState();
    double normalizer=eCF.getNormalizer();
    ofstream outputC;
    outputC.open("coeff.txt");
    for(size_t s=0; s<states.size(); ++s){
      outputC << verbatimCast(states.getDet(s)) << " " << sqrt(norm(states.getCoeff(s)))/sqrt(normalizer) << endl; 
    }
    outputC.close();
    myfile1 << count << " " << energy << " " <<   " " << aveEnergy<< endl;
    int allowedNumChangeSign = int(basis.getSize()*0.1);
    cout << "percentage of allowed sign change= " <<  allowedNumChangeSign<< endl;
    //cout << "sign= " << sign << endl;
    //cout << "lastSign= " << lastSign << endl;
    cout << "number of sign changes= " << abs(sign-lastSign) << endl;
    cout << "Ave energy= " << aveEnergy<< endl;
    cout << "energy= " << energy<< endl;
    //cout << "Exact energy= " << eMin<< endl;
    cout << "list size= " << list.size()<< endl;
    //cout << "abs(energy - lastEnergy)= " << abs(energy - lastEnergy) << endl;
    //cout << "fabs(energy - lastEnergy)= " << fabs(energy - lastEnergy) << endl;
    //cout << "abs(lastAveEnergy - aveEnergy)= " << abs(energy - lastEnergy) << endl;
    //while (abs(lastAveEnergy - aveEnergy) < 0.0005 || abs(lastEnergy-energy) < 0.0001){
    //while (abs(energy - lastEnergy) < 0.0001){
    //trainRate*=1.001;
    if (abs(sign-lastSign) > allowedNumChangeSign){
      trainRate*=0.90;
      cout << "trainRate=" << trainRate << endl;
    }
    //if (abs(lastAveEnergy - aveEnergy) < 0.001){
    //  trainRate+=0.002;
    //  cout << "trainRatex2=" << trainRate << endl;
    //}
    }
  }
  myfile1.close();
}
