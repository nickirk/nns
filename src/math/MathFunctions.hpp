/*
 * MathFunctions.hpp
 * Created on 08.3.2018
 * Author: Ke Liao 
 */
#ifndef MathFunctions_DEFINED
#define MathFunctions_DEFINED

#include <math.h>
#include <cmath>
#include <random>

double NormalDistribution(double dummy);
double TanhPrime(double in);
double Tanh(double in); 
double Linear(double in);
double LinearPrime(double in);
double Rectifier(double in);
double RectifierPrime(double in);
double Arcsinh(double in);
double ArcsinhPrime(double in);
double Gaussian(double in);
double GaussianPrime(double in);
double GaussianAntiSym(double in);
double GaussianAntiSymPrime(double in);
double Sigmoid(double in);
double SigmoidPrime(double in);
#endif
