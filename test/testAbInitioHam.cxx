/*
 * testAbInitioHam.cxx
 *
 *  Created on: Feb 13, 2018
 *      Author: Lauretta Schwarz
 */
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <numeric>
#include <math.h>
#include <stdio.h>
#include "../src/Nnw.hpp"
#include "../src/Basis.hpp"
#include "../src/Determinant.hpp"
#include "../src/AbInitioHamiltonian.hpp"
#include "../src/FermionicHamiltonian.hpp"
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
  int nexcit,nsingleexcit,ndoubleexcit;
  
  vector<int> spinConfig{spinUp, spinDown, numStates};
  int numHidden(10*numSites);
  //int numHidden1(2*numSites);
  vector<int> size_NNW = {numStates, numHidden, 2};
  //cout << "input number of hidden neurons=";
  //cin >> numHidden;
  bool readFromFile{false};
  double trainRate(1.5);
  //cout << "input training rate=";
  //cin >> trainRate;
  //generate basis, the basis class constructor takes in the spin configurations.
  Basis basis(spinConfig);
  //generate hamiltonian
  AbInitioHamiltonian modelHam(numStates);
  //double U{2.}, t{-1};
  string file_name = "../run/FCIDUMP"; 
  modelHam = readAbInitioHamiltonian(numStates, file_name);
  //// test Hubbard Hamiltonian
  //FermionicHamiltonian modelHam(numStates);
  //double U{8.0}, t{1.0};
  //modelHam = generateFermiHubbard(numStates, U, t);
  
  cout << "Basis size= " << basis.getSize() << endl;
  cout << "Determinants: " << endl;
  for (int i=0; i<basis.getSize(); ++i){
      std::vector<int> positions = getOccupiedPositions(basis.getDetByIndex(i));
      cout << " ";
      for (int j=0; j<positions.size(); ++j){
          cout << (positions[j]+1) << " , ";
      }
      cout << endl;
  }

  // form Hamiltonian and diagonalise it
  // Eigen matrices are columnmajor by default
  Eigen::MatrixXd H(Eigen::MatrixXd::Zero(basis.getSize(),
          basis.getSize()));
  for (int i=0; i<basis.getSize(); ++i){
      for (int j=0; j<basis.getSize(); ++j){
          double value = 0.0;
          value = modelHam(basis.getDetByIndex(j),basis.getDetByIndex(i));
          H(j,i) = value;
      }
  }
  cout << "\n Non-zero Hamiltonian matrix elements: i,j,H_ij \n" << endl;

  cout.precision(12);
  for (int i=0; i<basis.getSize(); ++i){
      for (int j=0; j<basis.getSize(); ++j){
          if (std::fabs(H(j,i))>1e-10){
              cout.width(10);
              cout << j; 
              cout.width(10); 
              cout << i;
              cout.width(25);
              std::fixed; 
              cout << H(j,i) << endl;
          }
      }
  }

  // write them to a file
  std::ofstream outfile;
  outfile.open("hamiltonian.txt",ios::out);
  outfile << "# i, j, H_ij" << endl;
  for (int j=0; j<basis.getSize(); ++j){
      for (int i=j; i<basis.getSize(); ++i){
          if (std::fabs(H(j,i))>1e-10){
              outfile.width(10);
              outfile << j+1; 
              outfile.width(10); 
              outfile << i+1;
              outfile.width(25);
              outfile.precision(12);
              std::fixed; 
              outfile << H(j,i) << endl;
          }
      }
  }
  outfile.close();

  // diagonalisation
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(H);
  if (eigensolver.info() != Success){
      abort();
  }
  cout << "\n The eigenvalues of H are:\n" << eigensolver.eigenvalues() << endl;
  //cout << "\n Matrix whose columns are eigenvector of H: \n" << eigensolver.eigenvectors() << endl;
  Eigen::MatrixXd eigenvecs(Eigen::MatrixXd::Zero(basis.getSize(),basis.getSize()));
  eigenvecs = eigensolver.eigenvectors();
  cout << "\n The eigenvectors of H:" << endl;
  for (int j=0; j<basis.getSize(); ++j){
      cout << "\n Eigenvector: " << (j+1) << "\n" << endl;
      for (int i=0; i< basis.getSize(); ++i){
          cout << i << "   " << eigenvecs(i,j) << endl;
      }
  }
  // write eigenvalues to a file
  Eigen::VectorXd eigenvals(Eigen::VectorXd::Zero(basis.getSize()));
  eigenvals = eigensolver.eigenvalues();
  outfile.open("eigenvalues.txt",ios::out);
  outfile << "# i, l_i" << endl;
  for (int j=0; j<eigenvals.size(); ++j){
      outfile.width(20);
      outfile << j+1; 
      outfile.width(23);
      outfile.precision(12);
      std::fixed; 
      outfile << eigenvals[j] << endl;
  }
  outfile.close();


  // test for excitation generators
  cout << "************************************************" << endl;
  cout << "\n Testing the excitation generators\n" << endl;
  cout << "************************************************" << endl;
  cout << "\n Deterministic excitation generator\n" << endl;
  cout << "************************************************" << endl;
  for (int i=0; i<basis.getSize(); ++i){
      cout << "\n \n Considering the source determinant:" << endl;
      std::vector<int> positions = getOccupiedPositions(basis.getDetByIndex(i));
      cout << "   ";
      for (int j=0; j<positions.size(); ++j){
          cout << (positions[j]+1) << " , ";
      }
      cout << endl;
      // count the number of excitations
      nexcit = 0;
      nsingleexcit = 0;
      ndoubleexcit = 0;
      nexcit = modelHam.countNumberCoupledStates(basis.getDetByIndex(i),3,nsingleexcit,ndoubleexcit);
      cout << "\n Number of single excitations: " << nsingleexcit << endl;
      cout << " Number of double excitations: " << ndoubleexcit << endl;
      cout << " Total number of excitations: " << nexcit << endl;
      std::vector<detType> CoupledDeterminants = modelHam.getCoupledStates(basis.getDetByIndex(i));
      cout << "\n List of excited determinants: " << endl;
      for (int j=0; j<CoupledDeterminants.size(); ++j){
          std::vector<int> positions = getOccupiedPositions(CoupledDeterminants[j]);
          cout << "   ";
          for (int k=0; k<positions.size(); ++k){
              cout << (positions[k]+1) << " , ";
          }
          cout << endl;
      }
      cout << "\n Predicted total number of excitations: " << nexcit << endl;
      cout << "\n Counted total number of excitations: " << CoupledDeterminants.size() << endl;
      if (nexcit != CoupledDeterminants.size()){
          cout << "\n \n !!!!!! Error! !!!!! \n" << endl;
          cout << " Predicted and Counted Number of excitations don't agree ! \n \n" << endl;
          cout << " Predicted number: " << nexcit << endl;
          cout << " Counted number: " << CoupledDeterminants.size() << endl;
      }
      int noffdiag = 0;
      for (int j=0; j<basis.getSize(); ++j){
          if (j==i) {continue;}
          if (std::fabs(H(j,i))>1e-10){
              noffdiag += 1;
          }
      }
      cout << "\n Counted total number of off-diagonal Hamiltonian matrix elements: " << noffdiag << endl;
      if ((nexcit < noffdiag)||(CoupledDeterminants.size() < noffdiag)){
          cout << "\n \n !!!!!! Error! !!!!! \n" << endl;
          cout << " Counted number of off-diagonal Hamiltonian matrix elements is larger than counted and / or predicted number of excitations ! \n \n" << endl;
          cout << " Predicted number: " << nexcit << endl;
          cout << " Counted number: " << CoupledDeterminants.size() << endl;
          cout << " Number of off-diagonal Hamiltonian matrix elements: " << noffdiag << endl;
      }
      if ((nexcit > noffdiag)||(CoupledDeterminants.size() > noffdiag)){
          cout << "\n \n !!!!!! Warning! !!!!! \n" << endl;
          cout << " Counted number of off-diagonal Hamiltonian matrix elements is smaller than counted and / or predicted number of excitations ! \n \n" << endl;
          cout << " Are some Hamiltonian matrix elements zero ?? If so, this is OK. \n \n" << endl;
          cout << " Predicted number: " << nexcit << endl;
          cout << " Counted number: " << CoupledDeterminants.size() << endl;
          cout << " Number of off-diagonal Hamiltonian matrix elements: " << noffdiag << endl;
      }
  }
  
  cout << "************************************************" << endl;
  cout << "\n Random excitation generator\n" << endl;
  cout << "************************************************" << endl;
  for (int i=0; i<basis.getSize(); ++i){
      cout << "\n \n Considering the source determinant:" << endl;
      cout << "   ";
      std::vector<int> positions = getOccupiedPositions(basis.getDetByIndex(i));
      for (int j=0; j<positions.size(); ++j){
          cout << (positions[j]+1) << " , ";
      }
      cout << endl;
      // count the number of excitations
      nexcit = 0;
      nsingleexcit = 0;
      ndoubleexcit = 0;
      nexcit = modelHam.countNumberCoupledStates(basis.getDetByIndex(i),3,nsingleexcit,ndoubleexcit);
      cout << "\n Number of single excitations: " << nsingleexcit << endl;
      cout << " Number of double excitations: " << ndoubleexcit << endl;
      cout << " Total number of excitations: " << nexcit << endl;
      std::vector<detType> CoupledDeterminants = modelHam.getCoupledStates(basis.getDetByIndex(i));
      std::sort(CoupledDeterminants.begin(),CoupledDeterminants.end(),compare_det);

      // histograms
      //std::vector<int> singlescount;
      //std::vector<double> singleshist;
      //std::vector<int> doublescount;
      //std::vector<double> doubleshist;
      std::vector<int> excitscount;
      std::vector<double> excitshist;

      //singlescount.assign(nsingleexcit, 0);
      //singleshist.assign(nsingleexcit, 0.0);
      //doublescount.assign(ndoubleexcit, 0);
      //doubleshist.assign(ndoubleexcit, 0.0);
      excitscount.assign(nexcit, 0);
      excitshist.assign(nexcit, 0.0);

      int counter = 1;
      double pgen=0.0;
      int iterations = 4000;
      while (counter <= iterations){

          pgen = 0.0;
          detType target = modelHam.getRandomCoupledState(basis.getDetByIndex(i),pgen);

          // find the determinant in the list and get its position
          int pos = std::find(CoupledDeterminants.begin(),CoupledDeterminants.end(),target) - 
              CoupledDeterminants.begin();

          if (pos < CoupledDeterminants.size()){
              // found the determinant
              excitscount[pos] += 1;
              excitshist[pos] += (1.0/pgen);

          }
          else{
              cout << "\n \n !!!!!! Error! !!!!! \n" << endl;
              cout << " Randomly generated determinant could not be found in list of connected determinants ! \n \n" << endl;
              cout << "\n \n Randomly generated determinant:" << endl;
              std::vector<int> positions = getOccupiedPositions(target);
              for (int j=0; j<positions.size(); ++j){
                  cout << (positions[j]+1) << " , ";
              }
              cout << endl;
          }

          // second way of evaluation generation probability
          double pgen_2 = modelHam.calcGenProp(basis.getDetByIndex(i),target);

          if (std::fabs(pgen-pgen_2) > 1e-10){
              cout << "\n \n !!!!!! Error! !!!!! \n" << endl;
              cout << " Evaluated generation probabilities don't agree ! \n \n" << endl;
              cout << "Probability 1: " << pgen << endl;
              cout << "Probability 2: " << pgen_2 << endl;
          }

          counter += 1;
      }

      int ndets = std::accumulate(excitscount.begin(), excitscount.end(), 0); 
      std::cout << "\n Number of randomly generated excitations: " << ndets << std::endl;

      std::string det_file = "hist_determinant_";
      det_file += std::to_string(i);
      outfile.open(det_file,ios::out);
      outfile << "# Number of target determinant, 1/pgen" << endl;
      double dev = 0.0;
      double devsign = 0.0;
      int nsampled = 0;
      for (int j=0; j<CoupledDeterminants.size(); ++j){
          double val = static_cast<double>(excitshist[j])/static_cast<double>(iterations);
          outfile.width(20);
          outfile << j+1; 
          outfile.width(23);
          outfile.precision(12);
          std::fixed; 
          outfile << val << endl;
          dev += std::fabs(val - 1.0);
          devsign += (val - 1.0);
          if (excitscount[j] != 0){
              nsampled += 1;
          }
      }
      outfile.close();

      std::cout << "\n Average absolute deviation from 1.0: " << 
          (dev/static_cast<double>(nsampled)) << std::endl;
      std::cout << "\n Average signed deviation from 1.0: " << 
          (devsign/static_cast<double>(nsampled)) << std::endl;

      if (nsampled != nexcit){
          std::cout << "\n !!! Warning !!! " << std::endl;
          std::cout << " Not all excitations accounted for... " << std::endl;
          std::cout << "Total number of excitations: " << nexcit << std::endl;
          std::cout << "Number of sampled excitations: " << nsampled << std::endl;
      }
  }
 
}