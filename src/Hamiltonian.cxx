/*
 * Hamiltonian.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */

#include <stdio.h>
#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>      // std::istringstream
#include <string>
#include <math.h>
#include <algorithm>
#include "Hamiltonian.hpp"

struct gathered{
	double val;
	int col,row;
};

Hamiltonian::Hamiltonian(double diagDelta_, double offDiagMagntd_, double
  offDiagNonZeroRatio_, Basis const &basis_): rowVec(std::vector<int>(0)),colVec(std::vector<int>(0)),
  valueVec(std::vector<double>(0)), basis(basis_){
      size = basis.getSize();
      diagDelta = diagDelta_;
      offDiagMagntd = offDiagMagntd_;
      offDiagNonZeroRatio = offDiagNonZeroRatio_;
      initHamiltonian();
      sortH();
      //std::cout <<"Value from within the class" <<std::endl; 
      //for (int i=0; i<valueVec.size(); ++i){
      //  std::cout << valueVec[i] << "," << std::endl;
      //}
      //for (int i=0; i<colVec.size(); ++i){
      //  std::cout << "Sorted column vec " << colVec[i] << "," << std::endl;
      //}
      //for (int i=0; i<rowVec.size(); ++i){
      //  std::cout << "Sorted row vec " << rowVec[i] << "," << std::endl;
      //}
}

int Hamiltonian::getSize() const {return size;}

int Hamiltonian::getSparseSize() const {return valueVec.size();}

void Hamiltonian::sparseAccess(int pos, int &row, int &col, double &value) const{
  row=rowVec[pos];col=colVec[pos];value=valueVec[pos];}

double Hamiltonian::operator() (detType const &i, detType const &j) const {
  double value(0.);
  int row = basis.getIndexByDet(i);
  int col = basis.getIndexByDet(j);
  value = (*this)(row, col);
  return value;
};

double Hamiltonian::operator() (int const i, int const j) const {
  std::vector<int>::const_iterator iter;
  double value(0.0);
  iter =  std::find(rowVec.begin()+lowerPos(j), rowVec.begin()+upperPos(j), i);
  int pos = iter - rowVec.begin();
  if (pos < upperPos(j)) value=valueVec[pos];
  //std::cout << "pos= " << pos << std::endl;
  //std::cout << "upperPos= " << upperPos(j)+1 << std::endl;
  return value;
}
//void Hamiltonian::sparseAcess(int pos, int &row, int &col, double &value){
//  row=rowVec[pos];col=colVec[pos];value=valueVec[pos];}

void Hamiltonian::multiplyMatrixVector(double *v, double *w){
    for(int i=0;i<size;++i){
      v[i]=0.0;
    }
    for(size_t i=0;i<valueVec.size();++i){
      v[rowVec[i]]+=valueVec[i]*w[colVec[i]];
    }
  }
      
void Hamiltonian::initHamiltonian(){
  //offDiagMagntd is the fraction of the diagDelta, it has value (0,1)
  //offDiagNonZeroRatio sets how much percentage of the offdiag terms are non0
  //diagonal terms
  for (int i = 0; i < size; i++){
    valueVec.push_back(0.5+diagDelta*pow(i,1.3));
    rowVec.push_back(i);
    colVec.push_back(i);
  }
  
  std::random_device rng;
  std::mt19937 eng(rng());
  std::uniform_real_distribution<> dist(0,1);
  //offdiagonal terms are generated sparsely and randomly
  for (int j = 0; j < size; j++){
    for (int i = j+1; i < size; i++){
      double prob = dist(eng);
      //std::cout << prob << std::endl;
      if (prob < offDiagNonZeroRatio){
        double fraction = dist(eng);
        double signProb = dist(eng);
        double sign(1);
        if (signProb < 0.5) sign = -1;
        //std::cout << fraction << std::endl;
        valueVec.push_back(sign*diagDelta*offDiagMagntd*fraction);
        valueVec.push_back(sign*diagDelta*offDiagMagntd*fraction);
        rowVec.push_back(i);
        rowVec.push_back(j);
        colVec.push_back(j);
        colVec.push_back(i);
      }
    }
  }
}

void Hamiltonian::sortH(){
	int numEls = getSparseSize();
	std::vector<gathered> buffer(numEls);
	for(int i=0;i<numEls;++i){
		buffer[i].val = valueVec[i];
		buffer[i].row = rowVec[i];
		buffer[i].col = colVec[i];
	}
	std::sort(buffer.begin(),buffer.end(),[]
			 (gathered const &a, gathered const &b){return a.col > b.col;});
	for(int i=0;i<numEls;++i){
		valueVec[i]=buffer[i].val;
		rowVec[i]=buffer[i].row;
		colVec[i]=buffer[i].col;
	}
}

int Hamiltonian::lowerPos(int const i) const{
  return std::distance(colVec.begin(),std::find(colVec.begin(),colVec.end(),i));
}

int Hamiltonian::upperPos(int const i) const{
  return std::distance(colVec.begin(),std::find(colVec.rbegin(),colVec.rend(),i).base());
}

void Hamiltonian::writeToFile(std::string Ham){
  std::ofstream ham;
  ham.open(Ham);
  for (int i=0; i<getSparseSize(); ++i){
    ham << valueVec[i] << " " << rowVec[i] << " " << colVec[i] <<std::endl;        
  } 
  ham.close();
  return;
}

void Hamiltonian::readFromFile(std::string Ham){
  rowVec.clear();
  colVec.clear();
  valueVec.clear();
  std::ifstream ham(Ham);
  std::string line;
  while (std::getline(ham, line)) {
    double value;
    int row, col;
    std::istringstream ss(line);
    ss >> value >> row >> col;
    valueVec.push_back(value);
    rowVec.push_back(row);
    colVec.push_back(col);
  }	
}
