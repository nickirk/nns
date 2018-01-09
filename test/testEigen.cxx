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
  Map<MatrixXd> m(pv+1,2,2);
  cout << "m=" << m << endl;
  v(2)=-1; 
  cout << "After m=" << m << endl;
}



