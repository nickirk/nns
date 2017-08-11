#include <stdio.h>
#include <vector>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <math.h>
using namespace Eigen;
class WaveBasis{
  public:
    int numEle;
    int numOrb;
    vector<vector<int>> occup;
    WaveBasis(): numEle(0), numOrb(0) {
    }
    WaveBasis(int numEle_, int numOrb_){
      numEle = numEle_;
      numOrb = numOrb_; 
    }
  private:
    void initWaveBasis(){
      
    }
   
}

inline int intcast(std::vector<bool> const & in){
  int out=0;
  for(size_t i=0;i<in.size();++i){
    if(in(i)){
      out+=pow(2,i);
    }
  }
}

class Hamiltonian{
public:					       
    int getSize() const {return size;}
    int getSparseSize() const {return valueVec.size();}
    void setDiagDelta(double newDiagDelta) {diagDelta = newDiagDelta;}
    void setOffDiagMagntd(double  newOffDiagMagntd) {offDiagMagntd = newOffDiagMagntd;}
    void setOffDiagNonZeroRatio(double newOffDiagNonZeroRatio) 
                       {offDiagNonZeroRatio = newOffDiagNonZeroRatio;}
    void setMatrixElements(double newDiagDelta, double newOffDiagMagntd, 
			 double newOffDiagNonZeroRatio) 
  {diagDelta = newDiagDelta; offDiagMagntd = newOffDiagMagntd; offDiagNonZeroRatio
  = newOffDiagNonZeroRatio;}
  void sparseAcess(int pos, int &row, int &col, double &value){
    row=rowVec[pos];col=colVec[pos];value=valueVec[pos];}
  double & operator() (std::vector<bool> i, std::vector<bool> j) {return entry(intcast(i),intcast(j));}
  double& operator() (int i, int j);
    Hamiltonian(int size_, double diagDelta_, double offDiagMagntd_, double
                offDiagNonZeroRatio_): rowVec(std::vector<int>(0)),colVec(std::vector<int>(0)),
                                       valueVec(std::vector<double>(0))
  {
      size = size_;
      diagDelta = diagDelta_;
      offDiagMagntd = offDiagMagntd_;
      offDiagNonZeroRatio = offDiagNonZeroRatio_;
      initHamiltonian();
    }
  void multMV(double *v, double *w){
    for(int i=0;i<size;++i){
      v[i]=0.0;
    }
    for(int i=0;i<valueVec.size();++i){
      v[rowVec[i]]+=valVec[i]*w[colVec[i]];
    }
  }
      
  private:
    int size;
    double diagDelta;
    double offDiagMagntd;
    double offDiagNonZeroRatio;
    std::vector<int> rowVec, colVec;
    std::vector<double> valueVec;
    void initHamiltonian(){
    //offDiagMagntd is the fraction of the diagDelta, it has value (0,1)
    //offDiagNonZeroRatio sets how much percentage of the offdiag terms are non0
      //diagonal terms
      for (int i = 0; i < size; i++){
        valueVec.push_back(0.0+diagDelta*i*i);
	rowVec.push_back(i);
	colVec.push_back(j);
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
};
