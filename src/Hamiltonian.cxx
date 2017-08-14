#include <stdio.h>
#include <vector>
#include <random>
#include <math.h>
#include "Hamiltonian.hpp"

struct gathered{
	double val;
	int col,row;
};

Hamiltonian::Hamiltonian(int size_, double diagDelta_, double offDiagMagntd_, double
  offDiagNonZeroRatio_): rowVec(std::vector<int>(0)),colVec(std::vector<int>(0)),
  valueVec(std::vector<double>(0)), sorted(false){
      size = size_;
      diagDelta = diagDelta_;
      offDiagMagntd = offDiagMagntd_;
      offDiagNonZeroRatio = offDiagNonZeroRatio_;
      initHamiltonian();
}

int Hamiltonian::getSize() const {return size;}

int Hamiltonian::getSparseSize() const {return valueVec.size();}

void Hamiltonian::sparseAccess(int pos, int &row, int &col, double &value){
  row=rowVec[pos];col=colVec[pos];value=valueVec[pos];}

double& Hamiltonian::operator() (detType i, detType j);

void Hamiltonian::multiplyMatrixVector(double *v, double *w){
    for(int i=0;i<size;++i){
      v[i]=0.0;
    }
    for(int i=0;i<valueVec.size();++i){
      v[rowVec[i]]+=valueVec[i]*w[colVec[i]];
    }
  }
      
void Hamiltonian::initHamiltonian(){
  //offDiagMagntd is the fraction of the diagDelta, it has value (0,1)
  //offDiagNonZeroRatio sets how much percentage of the offdiag terms are non0
  //diagonal terms
  for (int i = 0; i < size; i++){
    valueVec.push_back(0.0+diagDelta*i*i);
    rowVec.push_back(i);
    colVec.push_back(i);
  }
  
  //offdiagonal terms are generated sparsely and randomly
  for (int j = 0; j < size; j++){
    for (int i = j+1; i < size; i++){
      double prob = ((double) rand() / (RAND_MAX));
      //std::cout << prob << std::endl;
      if (prob < offDiagNonZeroRatio){
        double fraction = ((double) rand() / (RAND_MAX));
        //std::cout << fraction << std::endl;
        valueVec.push_back(diagDelta*offDiagMagntd*fraction);
        valueVec.push_back(diagDelta*offDiagMagntd*fraction);
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
	for(size_t i=0;i<numEls;++i){
		buffer[i].val = valueVec[i];
		buffer[i].row = rowVec[i];
		buffer[i].col = colVec[i];
	}
	std::sort(buffer.begin(),buffer.end(),[]
			 (gathered const &a, gathered const &b){return a.row > b.row;});
	for(size_t i=0;i<numEls;++i){
		valueVec[i]=buffer[i].val;
		rowVec[i]=buffer[i].row;
		colVec[i]=buffer[i].col;
	}
	sorted = true;
}
