#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <typeinfo>
using namespace Eigen;
using namespace std;
int main(){
  VectorXd v(9);
  v << 1,2,3,4,5,6,7,8,9;
  double *pv = &v(0);
  vector<Map<MatrixXd>> M;
  for(int i(0); i<2; ++i){
    Map<MatrixXd> m(pv+i,2,2);
    M.push_back(m);
  }
  M[0](0,0)=-1;
  M[1](0,0)=-3;
  cout << "v=" << endl;
  cout << v << endl;
  v(1)=-2;
  v(3)=-100;
  cout << M[0] << endl;
  cout << M[1] << endl;
  v(4) +=100;
  cout << M[0] << endl;
  cout << M[1] << endl;

  VectorXcd x(2);
  x(0)= std::complex<double>(1.,0);
  x(1)= std::complex<double>(0.,1.);
  MatrixXcd X;
  X = (x*x.adjoint()).adjoint();
  std::cout << "X=" << std::endl;
  std::cout << X << std::endl;
}




