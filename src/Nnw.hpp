#ifndef NeuralNetwork_DEFINED
#define NeuralNetwork_DEFINED

#include <vector>
#include <stdio.h>
#include <math.h>
#include <Eigen/Dense>
#include "Determinant.hpp"
#include "Basis.hpp"
#include "Hamiltonian.hpp"
#include "Sampler.hpp"
using namespace Eigen;
class NeuralNetwork{
  public:
    NeuralNetwork();
    NeuralNetwork(std::vector<int> sizes_);
    void train();

  private:
    std::vector<int> sizes;
    std::vector<VectorXd> activations;
    std::vector<VectorXd> biases;
    std::vector<MatrixXd> weights;
    std::vector<double> output_Cs;
    std::vector<VectorXd> delta_biases;
    std::vector<MatrixXd> delta_weights;
    double calcEnergy();
    double feedForward(double activation_);
    double NablaE_C();
    void backPropagate();
};

inline double tanh_prime(double in){return 1-tanh(in)*tanh(in);}
#endif
