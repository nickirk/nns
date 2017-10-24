#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include "../src/Nnw.hpp"
#include "../src/Basis.hpp"
#include "../src/Determinant.hpp"
#include "../src/Hamiltonian.hpp"
#include "../src/Sampler.hpp"
using namespace Eigen;

using namespace std;
int main(){
  int numSites(6);
  cout << "input number of sites=";
  cin >> numSites;
  int numStates(2*numSites);
  int numEle(4);
  cout << "input number of Ele=";
  cin >> numEle;
  int numHidden(15);
  int numHidden1(3);
  vector<int> size_NNW = {numStates, numHidden, 2};
  //cout << "input number of hidden neurons=";
  //cin >> numHidden;
  bool readFromFile{false};
  cout << "input readFromFile=";
  cin >> readFromFile;
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

  Eigen::MatrixXd H(Eigen::MatrixXd::Zero(basis.getSize(),basis.getSize()));
  for (int i(0); i<basis.getSize(); ++i){
	  for(int j(0); j<basis.getSize(); ++j){
	  H(i,j) = modelHam(basis.getDetByIndex(i),basis.getDetByIndex(j));
          cout << H(i,j) << " " ;
//"H( "<< i << "," << j << ")=" <<
	  }
     cout << endl;
    }
  Eigen::EigenSolver<MatrixXd> eSolver(H);
  VectorXd eVals = eSolver.eigenvalues().real();
  double eMin = 0.0;
  int pos = 0;
  for(int i = 0; i < basis.getSize(); ++i){
	  if((eVals(i)-eMin)<1e-8){
              cout << eVals(i) << " " << i << endl;
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
    myfilevec << eVector << endl;
  }
  vector<detType> list;
  for (int i=0; i< basis.getSize(); ++i){
    list.push_back(basis.getDetByIndex(i));
    vector<int> pos=getOccupiedPositions(basis.getDetByIndex(i));
    for (size_t j=0; j<pos.size(); j++){
      cout << pos[j] << ",";
    }
    cout << endl;
  }

  std::cout<<"Listsize= "<<list.size()<<std::endl;
  for(int i = 0; i< list.size(); ++i){
    std::cout<<"intCast= "<<verbatimCast(list[i])<<std::endl;
  }
  NeuralNetwork NNW(size_NNW, modelHam, basis);
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
  sampler.generateList(list); 
  sampler.removeDuplicate(list);
  //for (size_t i=0; i<list.size(); ++i){
  //    cout<<"intCast= " << verbatimCast(list[i])<<endl;
  //  }
  while (true){
    //list = NNW.train(list, 0.1); 
    //cout << "seeds size= " << list.size() << endl;
    //for (size_t i=0; i<list.size(); ++i){
    //  cout<<"intCast= " << verbatimCast(list[i])<<endl;
    //}
    lastSign = sign;
    list=NNW.train(list, trainRate);
    //NNW.train(list, trainRate);
    sampler.setReference(list);
    //cout << "Ref list size= " << list.size()<< endl;
    //for (size_t i=0; i<list.size(); ++i){
    //  cout<<"Ref intCast= " << verbatimCast(list[i])<<endl;
    //}
    sampler.generateList(list); 
    sampler.removeDuplicate(list);
    cout << "New list size= " << list.size()<< endl;
    //for (size_t i=0; i<list.size(); ++i){
    //  cout<<"New intCast= " << verbatimCast(list[i])<<endl;
    //}
    //energy = NNW.getEnergy();
    //sign = NNW.getSign();
    count++;
    if (count1 < 50){
      totalenergy+=energy; 
      count1++;
      aveEnergy = totalenergy/double(count1);
    }
    else{
      totalenergy+=energy; 
      count1++;
      aveEnergy = totalenergy/double(count1);
      totalenergy=0.;
      //count1=0;
      //list=NNW.train(list, trainRate);
      //NNW.train(list, trainRate);
      //sampler.setReference(list);
      //sampler.generateList(list); 
      //sampler.removeDuplicate(list);
      //cout << "size of seeds= " << list.size() << endl;
      //for (size_t i=0; i<list.size(); ++i){
      //  cout<<"Seeds intCast= " << verbatimCast(list[i])<<endl;
      //}
      //if (list.size() > int(1.2*numDetsToTrain_)){ 
      // numDetsToTrain_+= int(0.1*numDetsToTrain_);//numDetsToTrain_;
      // cout << "numDets= " << numDetsToTrain_<<endl;
      // sampler.setNumStates(numDetsToTrain_);
      // cout << "setNumStates= " << sampler.getNumStates() << endl;
      //}
    }
    if(count%1 == 0){
    //cout << "sign = " << sign<< endl;
    myfile1 << count << " " << energy << " " <<   " " << aveEnergy<< endl;
    coeffs = NNW.getCs();
    std::vector<Eigen::VectorXd> nablaE = NNW.getEnergyDerivative(list);
    for(size_t i = 0; i < coeffs.size(); ++i){
      cout << "Coeff on " << i << "\t" << "(" << pow(coeffs[i][0],2)+pow(coeffs[i][1],2)<<")"
     << ", exact " << eVector(i) << ", dE " << nablaE[i][0] << ", " << nablaE[i][1]<<endl;
    }
    int allowedNumChangeSign = int(basis.getSize()*0.1);
    cout << "percentage of allowed sign change= " <<  allowedNumChangeSign<< endl;
    //cout << "sign= " << sign << endl;
    //cout << "lastSign= " << lastSign << endl;
    cout << "number of sign changes= " << abs(sign-lastSign) << endl;
    cout << "Ave energy= " << aveEnergy<< endl;
    cout << "energy= " << energy<< endl;
    cout << "Exact energy= " << eMin<< endl;
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
    while (false){
      list = NNW.train(list, trainRate);
      sampler.setReference(list);
      cout << "size of seeds= " << list.size() << endl;
      for (size_t i=0; i<list.size(); ++i){
        cout<<"Seeds intCast= " << verbatimCast(list[i])<<endl;
      }
      if (list.size() == numDetsToTrain_){ 
       numDetsToTrain_+= 1;//numDetsToTrain_;
       cout << "numDets= " << numDetsToTrain_<<endl;
       sampler.setNumStates(numDetsToTrain_);
       cout << "setNumStates= " << sampler.getNumStates() << endl;
      }
      energy = NNW.getEnergy();
      lastSign = sign;
      sign = NNW.getSign();
      sampler.generateList(list);
      cout << "sign = " << sign<< endl;
      //totalenergy+=energy; 
      count++;
      aveEnergy = energy/double(count);
    }
  }
  myfile1.close();
}
