/*
 * Nnw.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */
#include <vector>
#include <iostream>
#include <math.h>
#include <cmath>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <complex>
#include "State.hpp"
#include "Nnw.hpp"
#include "errors.hpp"
#include "NormCF.hpp"
#include "Solver.hpp"
//using namespace Eigen;
NeuralNetwork::NeuralNetwork(std::vector<int> const &sizes_, 
CostFunction const &externalCF):sizes(sizes_), cf(&externalCF), sl(Solver(0.5)){
  momentumDamping = 0.6;
  momentum = false;
  //initial para for Nesterov's accelerated gradient descent
  lamdaS1 = 0.;
  lamdaS =0.;
  gammaS = 0.;
  gammaS1 = 0.;
  learningRate = 0.;
  iteration = 0;
  int numLayersBiasesWeights = sizes.size()-1;
  //calculate the size of the array:
  numNNP=0;
  for (int layer(0); layer<numLayersBiasesWeights; ++layer){
    numNNP+=sizes[layer]*sizes[layer+1]+sizes[layer+1];
  }
  yS = Eigen::VectorXd::Zero(numNNP);
  yS1 = Eigen::VectorXd::Zero(numNNP);
  Egz2 = Eigen::VectorXd::Zero(numNNP);
  NNP = Eigen::VectorXd::Ones(numNNP);
  generlisedForcePrev = Eigen::VectorXd::Zero(numNNP);
  //get the address of NNP
  adNNP = &NNP(0);
  nablaNNP = Eigen::VectorXd::Zero(numNNP);
  //get the address of nablaNNP
  adNablaNNP = &nablaNNP(0);
  //map parameters to matrices and initalize them with normal distribution.
  //map nablaPara to matrices and vectors with 0.
  int startPoint=0;
  for (int layer(0); layer<numLayersBiasesWeights; ++layer){
    //Pay special attention to weights, it has sizes.size()-1 layers,
    //instead of sizes.size() layers. Especially when reference to which
    //layer, should be careful.
    //0th layer for biases and weights starting from the 1st layer of neurons;
    //e.g biases[0] represents the biases of neurons on the 1st layer and
    // weights[0] represents the weights of connections between the 0th and 1st
    // layers of neurons.
    Eigen::Map<Eigen::MatrixXd> weightsTmp(adNNP+startPoint,sizes[layer+1],
		                           sizes[layer]);
    //NNP is modified here.
    //weightsTmp /= weightsTmp.size();
    weightsTmp = weightsTmp.unaryExpr(&NormalDistribution);
    weights.push_back(weightsTmp);
    nablaWeights.push_back(
      Eigen::Map<Eigen::MatrixXd>(adNablaNNP+startPoint, sizes[layer+1],
      sizes[layer])
    );
    startPoint+=sizes[layer]*sizes[layer+1];
    Eigen::Map<Eigen::VectorXd> biaseTmp(adNNP+startPoint,sizes[layer+1]); 
    //biaseTmp /= biaseTmp.size();
    biaseTmp = biaseTmp.unaryExpr(&NormalDistribution);
    biases.push_back(biaseTmp);
    nablaBiases.push_back(
      Eigen::Map<Eigen::VectorXd>(adNablaNNP+startPoint, sizes[layer+1])
    );
    startPoint+=sizes[layer+1];
  }

  //initial activity signals. 
  activations.push_back(Eigen::VectorXd::Zero(sizes[0]));
  inputSignal.push_back(Eigen::VectorXd::Zero(sizes[0]));
  //listDetsToTrain.push_back(fullBasis.getDetByIndex(0));
  for (int layer=0; layer < numLayersBiasesWeights; ++layer){
    activations.push_back(Eigen::VectorXd::Zero(sizes[layer+1]));
    inputSignal.push_back(Eigen::VectorXd::Zero(sizes[layer+1]));
   /*
    gFactorBiases.push_back(Eigen::VectorXd::Ones(sizes[layer+1]));
    gFactorBiasesPrev.push_back(Eigen::VectorXd::Ones(sizes[layer+1]));

    gFactorWeights.push_back(Eigen::MatrixXd::Ones(sizes[layer+1], sizes[layer]));
    gFactorWeightsPrev.push_back(Eigen::MatrixXd::Ones(sizes[layer+1], sizes[layer]));
   */
  }
  //nablaWeightsPrev = nablaWeights;
  //nablaBiasesPrev = nablaBiases;

  // Assign the output state and the evaluator to 0
  outputState = State();
}

void NeuralNetwork::train(std::vector<detType> const &listDetsToTrain, double eta, int iteration_){
// The coefficients are stored in scope of the train method and then stored into the state
  //nablaWeightsPrev = nablaWeights;
  //nablaBiasesPrev = nablaBiases;
  learningRate = eta;
  iteration  = iteration_;
  int numDets(listDetsToTrain.size());
  std::vector<detType> listDetsToTrainFiltered;
  std::vector<coeffType > outputCs;
  std::vector<std::vector<coeffType>> coupledOutputCsEpochs;
  std::vector<coeffType> coupledOutputCs;
  std::vector<std::vector<detType>> coupledDetsEpochs;
  std::vector<detType> coupledDets;
  int numLayersNeuron(sizes.size());
  int numLayersBiasesWeights(sizes.size()-1);
  //reset nablaWeights and nablaBiases to 0 values;
  nablaNNP *= 0.;
  //clear input signals in each training.
  inputSignalEpochs.clear();
  activationsEpochs.clear();

  sl.setGamma(eta);
  double prandom{0.0};
  std::random_device rng;
  double const normalizer = static_cast<double>(rng.max());
  double prob{1.0};
  //feed the first det to NNW and get the coeff. 
  coupledDets = getCoupledStates(listDetsToTrain[0]);
  coupledDetsEpochs.push_back(coupledDets);
  coeffType coupledCoeff;
  for (size_t i=0; i<coupledDets.size(); ++i){
    feedForward(coupledDets[i]);
    coupledCoeff=coeffType(activations[numLayersNeuron-1][0],
                           activations[numLayersNeuron-1][1]);
    coupledOutputCs.push_back(coupledCoeff);
  }
  coupledOutputCsEpochs.push_back(coupledOutputCs);
  //activations and inputSignal have been updated.
  feedForward(listDetsToTrain[0]);
  coeffType coeff(activations[numLayersNeuron-1][0], activations[numLayersNeuron-1][1]);
  coeffType coeffPrev =coeff;
  inputSignalEpochs.push_back(inputSignal);
  activationsEpochs.push_back(activations);
  outputCs.push_back(coeff);
  listDetsToTrainFiltered.push_back(listDetsToTrain[0]);
  for (int epoch=1; epoch < numDets; ++epoch){
    //initial input layer
    coupledOutputCs.clear();
    coupledDets.clear();
    coupledDets = getCoupledStates(listDetsToTrain[epoch]);
    for (size_t i=0; i<coupledDets.size(); ++i){
      feedForward(coupledDets[i]);
      coupledCoeff=coeffType(activations[numLayersNeuron-1][0],
                             activations[numLayersNeuron-1][1]);
      coupledOutputCs.push_back(coupledCoeff);
    }
    feedForward(listDetsToTrain[epoch]);
    coeff = coeffType(activations[numLayersNeuron-1][0], activations[numLayersNeuron-1][1]); 
    prob=norm(coeff)/norm(coeffPrev);
    prandom=rng()/normalizer;
    coeffPrev = coeff;
    //if ((prob-1)>1.e-8 || (prandom-prob) < -1.e-8){
    if (true){
      inputSignalEpochs.push_back(inputSignal);
      activationsEpochs.push_back(activations);
      outputCs.push_back(coeff);
      coupledDetsEpochs.push_back(coupledDets);
      coupledOutputCsEpochs.push_back(coupledOutputCs);
      listDetsToTrainFiltered.push_back(listDetsToTrain[epoch]);
    } 
  }
  
  // State generation fails if number of determinants and coefficients are not equal
  
  try{
    outputState = State(listDetsToTrainFiltered, outputCs, coupledDetsEpochs, coupledOutputCsEpochs);
  }
  catch(sizeMismatchError const &ex){
          // Check if outputCs and coupledOutputCsEpochs have the same size.
	  // If the two are not equal in size, shrink them
	  int mSize = (listDetsToTrainFiltered.size() > outputCs.size())?outputCs.size():listDetsToTrain.size();
	  outputCs.resize(mSize);
	  listDetsToTrainFiltered.resize(mSize);
	  outputState = State(listDetsToTrainFiltered,outputCs, coupledDetsEpochs,
                              coupledOutputCsEpochs);
  }

  updateParameters(2);

}

void NeuralNetwork::updateParameters(int method){
  //method corresponds to
  // 0: Stochastic gradiend desend
  // 1: Stochastic reconfiguration
  // 2: Nesterov's Accelerated Gradient Descent
  Eigen::VectorXd generlisedForce=backPropagate(inputSignalEpochs, activationsEpochs);
  if (method == 0){
    sl.update(NNP,generlisedForce);
    //update weights and biases
  }
  //Stochastic reconfiguration
  else if (method == 1){
    Eigen::MatrixXcd dCdw=backPropagateSR(inputSignalEpochs, activationsEpochs);
    std::vector<coeffType> outputCs = outputState.getAllCoeff();
    Eigen::Map<Eigen::VectorXcd> ci(&(outputCs[0]),outputCs.size());
    sl.update(NNP,generlisedForce,ci,dCdw, iteration);
  }
  else if (method ==2){
     lamdaS1 = (1+std::sqrt(1+4*lamdaS*lamdaS))/2.;
     gammaS = (1-lamdaS)/(lamdaS1);//*std::exp(-1./100*iteration);
     //std::cout << "here" << std::endl;
     //std::cout << "generlisedForce" << std::endl;
     //std::cout << generlisedForce << std::endl;
     //std::cout << "gammaS="<< gammaS << std::endl;
     double rho=0.9;//*std::exp(-1./100*iteration);
     Egz2 = rho*Egz2.matrix() + (1-rho)*generlisedForce.array().square().matrix();
     Egz2 += Eigen::VectorXd::Ones(numNNP)*1e-4;
     Eigen::VectorXd RMS = (Egz2).array().sqrt();
     Eigen::VectorXd tau = learningRate * RMS.array().inverse();
     //yS1 = NNP - learningRate * generlisedForce;
     yS1 = (NNP.array() -  tau.array() * generlisedForce.array()).matrix();
     if (iteration == 0) yS=yS1;
     //if (iteration == 300) gammaS = 0.5;
     NNP = (1-gammaS)*yS1 + gammaS*yS;
     lamdaS = lamdaS1;
     yS = yS1; 
     std::cout << "gammaS" << std::endl;
     std::cout << gammaS << std::endl;
     std::cout << "lamdaS" << std::endl;
     std::cout << lamdaS << std::endl;
     std::cout << "tau" << std::endl;
     std::cout << tau << std::endl;
  }
  else if (method == 3){
    if (iteration==0) generlisedForcePrev = generlisedForce;
    std::cout << generlisedForce-generlisedForcePrev << std::endl;
    std::cout << "generlisedForcePrev" << iteration << std::endl;
    double beta = (generlisedForce.transpose()*(generlisedForce-generlisedForcePrev)/
                  (generlisedForcePrev.transpose()*generlisedForcePrev))(0);
    std::cout << (generlisedForce.transpose()*(generlisedForce-generlisedForcePrev)/
           (generlisedForcePrev.transpose()*generlisedForcePrev)) <<std::endl;
    std::cout << "beta=" << beta << std::endl;
    generlisedForce += beta*generlisedForcePrev; 
    Eigen::MatrixXcd dCdw=backPropagateSR(inputSignalEpochs, activationsEpochs);
    std::vector<coeffType> outputCs = outputState.getAllCoeff();
    Eigen::Map<Eigen::VectorXcd> ci(&(outputCs[0]),outputCs.size());
    //sl.update(NNP,generlisedForce,ci,dCdw, iteration);
    NNP -= learningRate * generlisedForce;
    generlisedForcePrev = generlisedForce;
  }
}

Eigen::VectorXd NeuralNetwork::feedForward(detType const& det) const{
  int numStates=det.size();
  int numLayersNeuron(sizes.size());
  for (int state=0; state<numStates; ++state){
    activations[0][state] = det[state]?1.0:0.0;
    inputSignal[0][state] = det[state]?1.0:0.0;
  }

  for (int layer=1; layer < numLayersNeuron; ++layer){
    //Here the 0th layer of weights correspond to the connections between
    //the 0th and 1st layer. Here layer refers to the Neuron layer. When use 
    //it to refer the Biases and weights' layer, we need conversion.
    activations[layer] = weights[layer-1]*activations[layer-1]+biases[layer-1];
    inputSignal[layer] = activations[layer];
    if (layer == numLayersNeuron-1)
    activations[layer] = activations[layer].unaryExpr(&Linear);
    else 
    activations[layer] = activations[layer].unaryExpr(&Tanh);
  }
  return activations[numLayersNeuron-1];
}


Eigen::VectorXd NeuralNetwork::backPropagate(
       std::vector<std::vector<Eigen::VectorXd>> const &inputSignalEpochs,
       std::vector<std::vector<Eigen::VectorXd>> const &activationsEpochs
     ){
  int numDets = outputState.size();
  int numLayersNeuron(sizes.size());
  int numLayersBiasesWeights(sizes.size()-1);
  //everytime the backPropagate is called, we should reset nabla* to zero.
  nablaNNP *= 0.;
  std::vector<Eigen::VectorXd> dEdC = cf->nabla(outputState);
  for (int epoch=0; epoch < numDets; ++epoch){
    Eigen::VectorXd deltaTheLastLayer;
    deltaTheLastLayer = 
    (dEdC[epoch].array() * 
    inputSignalEpochs[epoch][numLayersNeuron-1].unaryExpr(&Linear_prime).array()).matrix();
    //adding up all nablaBiases and nablaWeights, in the end use the average
    //of them to determine the final change of the weights and biases.
    //Treat them as Vectors and Matrices.
    //Starting from the last layer of neuron and backpropagate the error.
    //nablaWeights have the same structure as weights, only numLayers-2 layers.
    std::vector<Eigen::VectorXd> deltaAllLayers;
    deltaAllLayers.push_back(deltaTheLastLayer);
    nablaBiases[numLayersBiasesWeights-1] += deltaTheLastLayer;
    nablaWeights[numLayersBiasesWeights-1] += 
      deltaTheLastLayer * activationsEpochs[epoch][numLayersNeuron-2].transpose();
    for (int layer=numLayersBiasesWeights-2; layer >= 0; --layer){
      //Calculating the error from the second last layer to
      //the 1st layer (0th layer is the input neurons, have no biases.)
      //Remember weights have only 0th -- numLayers-2 layers.
      Eigen::VectorXd deltaPreviousLayer=deltaAllLayers[numLayersBiasesWeights-2-layer];
      Eigen::VectorXd deltaThisLayer;
      //MatrixXd deltaAsDiagonal;
      //\delta_l = ((weights_{l+1, l}^T\delta^{l+1})) .* tanh'(z^l);
      //where ^T means transpose, l means lth layer, tanh' means the derivative
      //of tanh and z^l= tanh(weights_{l,l-1}activations^{l-1}) is the input signal 
      //on the lth layer neurons, where weights_{l, l-1} corresponds to the 
      //connections between the lth and l-1 th layers of neurons. Notice
      // .* means elementwise multiply.
      
      deltaThisLayer = weights[layer+1].transpose() * deltaPreviousLayer; 
      //To achive elementwise multiply, form a diagonal vector.
      //deltaAsDiagonal = deltaThisLayer.asDiagonal();

      deltaThisLayer = deltaThisLayer.array()
         * inputSignalEpochs[epoch][layer+1].unaryExpr(&Tanh_prime).array();
      deltaAllLayers.push_back(deltaThisLayer);
      nablaBiases[layer] += deltaThisLayer;
      //get a weight matrix. \partial C/\partial w^{l}_{jk} = a^{l-1}_k \delta_j^l 
      //the layer here refers to the lth layer of Biases and weights, so for
      //activation layer refers to the l-1th layer.
      nablaWeights[layer] += deltaThisLayer * activationsEpochs[epoch][layer].transpose();
    }
  }
  return nablaNNP;
}

Eigen::MatrixXcd NeuralNetwork::backPropagateSR( std::vector<std::vector<Eigen::VectorXd>> const &inputSignalEpochs,
       std::vector<std::vector<Eigen::VectorXd>> const &activationsEpochs){
//This step produce a complex matrix dCdw. It is done via the same backPropagate
//algorithm, with dEdC set to a vector of (1,0) (real part) and a vector of (0,1)
//(imaginary part). 
  int numDets = outputState.size();
  int numLayersNeuron(sizes.size());
  int numLayersBiasesWeights(sizes.size()-1);
  //Here, not like in backPropagate,
  //the dEdC is just a parameter which we choose to produce the real
  //or imaginary part of dCdw, doesn't mean dEdC itself.
  Eigen::MatrixXcd dCdw(numNNP,numDets);
  Eigen::Vector2d dedc;
  Eigen::Vector2d mask;
  mask << -1., 1.;
  for (int i(0); i <=1; ++i){
  //loop through dedc = (1,0) and (0,1) which will give us the real and imaginary part of
  //dCdw, respectivly.
    Eigen::MatrixXd dCdwTmp(numNNP,numDets);
    dedc << 1, 0;
    dedc += mask * i;
    std::vector<Eigen::VectorXd> dEdC(numDets,dedc);
    //create w_k|--C_i matrix
    for (int epoch(0); epoch < numDets; ++epoch){
      //reset nablaNNP to 0
      nablaNNP *= 0.;
      Eigen::VectorXd deltaTheLastLayer;
      deltaTheLastLayer = 
      (dEdC[epoch].array() * 
      inputSignalEpochs[epoch][numLayersNeuron-1].unaryExpr(&Linear_prime).array()).matrix();
      //adding up all nablaBiases and nablaWeights, in the end use the average
      //of them to determine the final change of the weights and biases.
      //Treat them as Vectors and Matrices.
      //Starting from the last layer of neuron and backpropagate the error.
      //nablaWeights have the same structure as weights, only numLayers-2 layers.
      std::vector<Eigen::VectorXd> deltaAllLayers;
      deltaAllLayers.push_back(deltaTheLastLayer);
      nablaBiases[numLayersBiasesWeights-1] = deltaTheLastLayer;
      nablaWeights[numLayersBiasesWeights-1] = 
        deltaTheLastLayer * activationsEpochs[epoch][numLayersNeuron-2].transpose();
      for (int layer=numLayersBiasesWeights-2; layer >= 0; --layer){
        //Calculating the error from the second last layer to
        //the 1st layer (0th layer is the input neurons, have no biases.)
        //Remember weights have only 0th -- numLayers-2 layers.
        Eigen::VectorXd deltaPreviousLayer=deltaAllLayers[numLayersBiasesWeights-2-layer];
        Eigen::VectorXd deltaThisLayer;
        //MatrixXd deltaAsDiagonal;
        //\delta_l = ((weights_{l+1, l}^T\delta^{l+1})) .* tanh'(z^l);
        //where ^T means transpose, l means lth layer, tanh' means the derivative
        //of tanh and z^l= tanh(weights_{l,l-1}activations^{l-1}) is the input signal 
        //on the lth layer neurons, where weights_{l, l-1} corresponds to the 
        //connections between the lth and l-1 th layers of neurons. Notice
        // .* means elementwise multiply.
        deltaThisLayer = weights[layer+1].transpose() * deltaPreviousLayer; 
        deltaThisLayer = deltaThisLayer.array()
           * inputSignalEpochs[epoch][layer+1].unaryExpr(&Tanh_prime).array();
        deltaAllLayers.push_back(deltaThisLayer);
        nablaBiases[layer] = deltaThisLayer;
        //get a weight matrix. \partial C/\partial w^{l}_{jk} = a^{l-1}_k \delta_j^l 
        //the layer here refers to the lth layer of Biases and weights, so for
        //activation layer refers to the l-1th layer.
        nablaWeights[layer] = deltaThisLayer * activationsEpochs[epoch][layer].transpose();
      }
      //fill up the dCdwTmp matrix column by column. Column index is C_i, row index is w_k
      dCdwTmp.col(epoch) << nablaNNP; 
    }
    if (i==0) dCdw.real()=dCdwTmp;
    if (i==1) dCdw.imag()=dCdwTmp;
  }
  return dCdw;
}

void preTrain(NeuralNetwork &network, State const &target, double trainRate, int iteration){
// Trains the network to represent some state target
// We first backup the current cost function
	CostFunction const *backupCF = network.getCostFunction();
// Then, set the cost function to the L2-distance to target
	NormCF stateDistance(target);
	network.setCostFunction(stateDistance);
// Set up an initial list of determinants to train
// Caveat: All determinants not present in the state do not matter
// i.e. their coefficients are treated as unknown
	std::vector<detType > list = target.getDets();
// Train the network
	int const maxTrainCount = 1000;
	for(int i = 0; i < maxTrainCount;++i){
		network.train(list, trainRate, iteration);
		std::cout<<"Distance " << network.getEnergy() << std::endl;
	}
	network.setCostFunction(*backupCF);
}

double NormalDistribution(double input)
{
  //input will be a constant, which is a normalization factor.
  static std::mt19937 rng;
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
