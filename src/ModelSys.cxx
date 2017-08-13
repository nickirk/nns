#include <ModelSys.hpp>
#include <stdio.h>
#include <vector>
#include <random>
#include <math.h>


Basis::Basis(int numEle_, int numOrb_){
  numEle = numEle_;
  numOrb = numOrb_; 
  generateBasis();
}

int Basis::getSize(){return size;}

void Basis::generateBasis(){
  basis = new BasisIndexRef[size];
  Determinant tmpDet(numOrb);
  Determinant lastDet(numOrb);
  for (int i=numOrb-1; i>=(numOrb-numEle); i--){
    lastDet.create(i);
  }
  for (int i=0; i<numEle; i++){
    tmpDet.create(i);
  }  
  int cursor(numEle-1);
  int end(numOrb-1);
  while (tmpDet.intCast() != lastDet.intCast()){
    if (cursor == end){
      cursor--;
      end--;
    }  
  }
}

//detType Basis::getIndex(detType det_) {
//  std::vector<int> position;
//  for (int i=0; i < numOrb; i++){
//    if (det_[i] == 1) position.push_back(i);
//  }
//  for (int j=0; j<numEle; j++){
//    
//  }
//}


Determinant Basis::getDetByIndex(int index_){
  return 
}

class Basis::BasisIndexRef(){
  public:
    BasisIndexRef(){
       
      index = 0;
    }
    BasisIndexRef(Determinant det_, int index_): det(det_), index(index_){}
    //static bool sortByDet (BasisIndexRef const& basisindexref1, 
    //  BasisIndexRef const& basisindexref2) {
    //  
    //}
    static bool sortByIndex (BasisIndexRef const& basisindexref1, 
      BasisIndexRef const& basisindexref2) {
      return basisindexref1.index < basisindexref2.index;
    }
  private:
    Determinant det(size);
    int index;
}

Hamiltonian::Hamiltonian(int size_, double diagDelta_, double offDiagMagntd_, double
  offDiagNonZeroRatio_): rowVec(std::vector<int>(0)),colVec(std::vector<int>(0)),
  valueVec(std::vector<double>(0)){
      size = size_;
      diagDelta = diagDelta_;
      offDiagMagntd = offDiagMagntd_;
      offDiagNonZeroRatio = offDiagNonZeroRatio_;
      initHamiltonian();
}

int Hamiltonian::getSize() const {return size;}

int Hamiltonian::getSparseSize() const {return valueVec.size();}

void Hamiltonian::sparseAcess(int pos, int &row, int &col, double &value){
  row=rowVec[pos];col=colVec[pos];value=valueVec[pos];}

double& Hamiltonian::operator() (std::vector<bool> i, std::vector<bool> j) {return entry(intcast(i),intcast(j));}

double& Hamiltonian::operator() (int i, int j) {};

void Hamiltonian::multiplyMatrixVector(double *v, double *w){
    for(int i=0;i<size;++i){
      v[i]=0.0;
    }
    for(int i=0;i<valueVec.size();++i){
      v[rowVec[i]]+=valVec[i]*w[colVec[i]];
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