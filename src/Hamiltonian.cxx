#include <stdio.h>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <math.h>
using namespace Eigen;
class Hamiltonian{
  public:
    int size;
    double diagDelta;
    double offDiagMagntd;
    double offDiagNonZeroRatio;
    MatrixXd entry;
    Hamiltonian(int size_, double diagDelta_, double offDiagMagntd_, double
                offDiagNonZeroRatio_){
      size = size_;
      diagDelta = diagDelta_;
      offDiagMagntd = offDiagMagntd_;
      offDiagNonZeroRatio = offDiagNonZeroRatio_;
      initHamiltonian();
    }
  private:
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
          std::cout << prob << std::endl;
          if (prob < offDiagNonZeroRatio){
            double fraction = ((double) rand() / (RAND_MAX));
            std::cout << fraction << std::endl;
            entry(j, i) = diagDelta*offDiagMagntd*fraction;
            entry(i, j) = entry(i, j);
          }
        }
      }
    }
};
