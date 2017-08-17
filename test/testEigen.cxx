#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
using namespace Eigen;
int main(){
  VectorXd m(VectorXd::Random(3));
  VectorXd n(VectorXd::Random(2));
  
  std::cout << m*n.transpose();
}


