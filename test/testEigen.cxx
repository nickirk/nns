#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <typeinfo>
using namespace Eigen;
using namespace std;
int main(){
MatrixXd ones = MatrixXd::Ones(3,3);
MatrixXd ones1 = MatrixXd::Ones(3,3);
EigenSolver<MatrixXd> es(ones);
EigenSolver<MatrixXd> es1(ones1);
cout << "The first eigenvector of the 3x3 matrix of ones is:"
     << endl << es1.eigenvectors().col(0) << endl;
cout << "type of es1= " << typeid(es1.eigenvectors().col(0)).name() << endl;
auto eVector1=es.eigenvectors().col(0);
VectorXcd eVector=es.eigenvectors().col(0);
//VectorXcd eVector=eVector1;
cout << "type of eVector= " << typeid(eVector).name() << endl;
cout << "type of eVector1= " << typeid(eVector1).name() << endl;
auto eigenValue=es.eigenvalues();

while (true){
cout << "eVector= " << eVector << endl;
cout << "eVector1= " << eVector1 << endl;
cout << "eigenValues= " << eigenValue << endl;
}
}



