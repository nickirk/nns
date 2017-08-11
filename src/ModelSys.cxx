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
    int getSize() const {return size_;}
    void setDiagDelta(double newDiagDelta) {diagDelta = newDiagDelta;}
    void setOffDiagMagntd(double  newOffDiagMagntd) {offDiagMagntd = newOffDiagMagntd;}
    void setOffDiagNonZeroRatio(double newOffDiagNonZeroRatio) 
                       {offDiagNonZeroRatio = newOffDiagNonZeroRatio;}
    void setMatrixElements(double newDiagDelta, double newOffDiagMagntd, 
			 double newOffDiagNonZeroRatio) 
  {diagDelta = newDiagDelta; offDiagMagntd = newOffDiagMagntd; offDiagNonZeroRatio
  = newOffDiagNonZeroRatio;}
  double & operator() (std::vector<bool> i, std::vector<bool> j) {return entry(intcast(i),intcast(j));}
    Hamiltonian(int size_, double diagDelta_, double offDiagMagntd_, double
                offDiagNonZeroRatio_){
      size = size_;
      diagDelta = diagDelta_;
      offDiagMagntd = offDiagMagntd_;
      offDiagNonZeroRatio = offDiagNonZeroRatio_;
      initHamiltonian();
    }
  private:
    int size;
    double diagDelta;
    double offDiagMagntd;
    double offDiagNonZeroRatio;
    MatrixXd entry;
    void initHamiltonian(){
    //offDiagMagntd is the fraction of the diagDelta, it has value (0,1)
    //offDiagNonZeroRatio sets how much percentage of the offdiag terms are non0
      entry = MatrixXd::Zero(size, size);
      //diagonal terms
      for (int i = 0; i < size; i++){
        entry(i,i) = 0.0+diagDelta*i*i;  
      }
  
      //offdiagonal terms are generated sparsely and randomly
      for (int j = 0; j < size; j++){
        for (int i = j+1; i < size; i++){
          double prob = ((double) rand() / (RAND_MAX));
          //std::cout << prob << std::endl;
          if (prob < offDiagNonZeroRatio){
            double fraction = ((double) rand() / (RAND_MAX));
            //std::cout << fraction << std::endl;
            entry(j, i) = diagDelta*offDiagMagntd*fraction;
            entry(i, j) = entry(j, i);
	             }
        }
      }
    }
};
