/*
 * MathFunctions.cxx
 * Created on 08.3.2018
 * Author: Ke Liao 
 */

double NormalDistribution(double input)
{
  //input will be a constant, which is a normalization factor.
  std::random_device rd;
  static std::mt19937 rng(rd());
  static std::normal_distribution<> nd(0.,input);
  return nd(rng);
}

double Tanh_prime(double in){return 1-tanh(in)*tanh(in);};
double Tanh(double in){return tanh(in);};
double Linear(double in) {return in;};
double Linear_prime(double in){return 1;};
double Rectifier(double in){return (in<-1e-8)?0.01*in:in;};
double Rectifier_prime(double in){return (in<-1e-8)?0.01:1.;};
double Arcsinh(double in){return std::asinh(in);};
double Arcsinh_prime(double in){return 1./sqrt(1.+ in*in);};
double Gaussian(double in){return exp(-pow(in,2));};
double Gaussian_prime(double in){return -2*in*Gaussian(in);};
double GaussianAntiSym(double in){
  double value(0.);
  if (in > 1e-8) value = 1.-Gaussian(in);
  else if (in < -1e-8) value = Gaussian(in)-1.;
  return value;
}
double GaussianAntiSym_prime(double in){
  double value(0.);
  if (in > 1e-8) value = -Gaussian_prime(in);
  else if (in < -1e-8) value = Gaussian_prime(in);
  return value;
}
double Sigmoid(double in){return 1./(1+exp(-in));};
double Sigmoid_prime(double in){return Sigmoid(in)*(1-Sigmoid(in));};
