/*
 * MathFunctions.hpp
 * Created on 08.3.2018
 * Author: Ke Liao 
 */
#ifndef MathFunctions_DEFINED
#define MathFunctions_DEFINED

#include <math.h>
#include <cmath>

double NormalDistribution(double dummy);
double Tanh_prime(double in);
double Tanh(double in); 
double Linear(double in);
double Linear_prime(double in);
double Rectifier(double in);
double Rectifier_prime(double in);
double Arcsinh(double in);
double Arcsinh_prime(double in);
double Gaussian(double in);
double Gaussian_prime(double in);
double GaussianAntiSym(double in);
double GaussianAntiSym_prime(double in);
double Sigmoid(double in);
double Sigmoid_prime(double in);

#endif
