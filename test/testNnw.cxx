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
  int numSites(3);
  int numStates(2*numSites);
  int numEle(3);
  int numHidden(15);
  vector<int> size_NNW = {numStates, numHidden, 2};
  bool readFromFile{false};
  double trainRate(0.005);
  Basis basis(numStates,numEle);
  Hamiltonian modelHam(numStates);
  double U{4}, t{-1};
  modelHam = generateHubbard(numStates, U, t);
  cout << "Basis size= " << basis.getSize() << endl;
  cout << "print out Ham element" << endl;

  // First, do an exact diagonalization to check for the exact energy
  // therefore, construct the Hamiltonian stored in a dense matrix format
  Eigen::MatrixXd H(Eigen::MatrixXd::Zero(basis.getSize(),basis.getSize()));
  for (int i(0); i<basis.getSize(); ++i){
	  for(int j(0); j<basis.getSize(); ++j){
	  H(i,j) = modelHam(basis.getDetByIndex(i),basis.getDetByIndex(j));
	  }
    }
  // now, get the eigenvalues
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

  // Write all eigenvalues to a file
  if(!readFromFile){
    ofstream myfile("eigen.txt");
    cout << "writing eigen values to file" << endl;
    myfile << eVals <<  endl;

    myfile.close();
  } 
  vector<detType> list;
  for (int i=0; i< basis.getSize(); ++i){
    list.push_back(basis.getDetByIndex(i));
    vector<int> pos=getOccupiedPositions(basis.getDetByIndex(i));
  }

  // Set up the Nnw
  std::cout<<"Listsize= "<<list.size()<<std::endl;
  NeuralNetwork NNW(size_NNW, modelHam, basis);
  detType HF=basis.getDetByIndex(0);
  ofstream myfile1;
  myfile1.open ("energy.txt");
  double aveEnergy(0.);
  double totalenergy(0.);
  double energy(0.);
  int sign(0);
  int lastSign(0);
  int count1(0);
  std::vector<Eigen::VectorXd> coeffs;
  int const maxCount = 7000;
  for(int count = 0; count < maxCount; ++count){
    lastSign = sign;
    NNW.train(list, trainRate);
    energy = NNW.getEnergy();
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
    }
    if(count%100 == 0){
    myfile1 << count << " " << energy << " " <<   " " << aveEnergy<< endl;
    int allowedNumChangeSign = int(basis.getSize()*0.1);
    cout << "percentage of allowed sign change= " <<  allowedNumChangeSign<< endl;
    cout << "number of sign changes= " << abs(sign-lastSign) << endl;
    cout << "Ave energy= " << aveEnergy<< endl;
    cout << "energy= " << energy<< endl;
    cout << "Exact energy= " << eMin<< endl;
    cout << "list size= " << list.size()<< endl;
    if (abs(sign-lastSign) > allowedNumChangeSign){
      trainRate*=0.90;
      cout << "trainRate=" << trainRate << endl;
    }
    }
  }
  myfile1.close();
}
