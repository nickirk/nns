#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <typeinfo>
using namespace Eigen;
using namespace std;
int main(){
VectorXd m = VectorXd::Random(3);

cout << "m[]=" << m(1) << endl;
}



