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
//using namespace Eigen;
NeuralNetwork::NeuralNetwork(std::vector<int> const &sizes_, CostFunction const &externalCF):sizes(sizes_), cf(&externalCF), sl(Solver(1)){
  //momentumDamping = 0.6;
  //momentum = false;
  //epsilon = 0.5;
  int numLayersBiasesWeights = sizes.size()-1;
  //calculate the size of the array:
  int numPara(0);
  for (int layer(0); layer<numLayersBiasesWeights; ++layer){
    numPara+=sizes[layer]*sizes[layer+1]+sizes[layer];
  }
  NNP = Eigen::VectorXd::Ones(NNP);
  //get the address of NNP
  adNNP = &NNP(0);
  nablaNNP = Eigen::VectorXd::Zero(NNP);
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
    Eigen::Map<Eigen::MatrixXd> weightsTmp(adNNP+startPoint,sizes[layer+1],sizes[layer])
    //NNP is modified here.
    weightsTmp /= weightsTmp.size();
    weightsTmp = weightsTmp.unaryExpr(&NormalDistribution);
    weights.push_back(weightsTmp);
    nablaWeights.push_back(
      Eigen::Map<Eigen::MatrixXd>(adNablaNNP+startPoint, sizes[layer+1],
      sizes[layer])
    );
    startPoint+=sizes[layer]*sizes[layer+1];
    Eigen::Map<Eigen::VectorXd> biaseTmp(adNNP+startPoint,sizes[layer+1]); 
    biaseTmp /= biaseTmp.size();
    biaseTmp = biaseTmp.unaryExpr(&NormalDistribution);
    biases.push_back(biaseTmp);
    nablaWeights.push_back(
      Eigen::Map<Eigen::VectorXd>(adNablaNNP+startPoint, sizes[layer+1]);
    );
    startPoint+=sizes[layer]*sizes[layer+1]+sizes[layer];
  }

  //initial activity signals. 
  activations.push_back(Eigen::VectorXd::Zero(sizes[0]));
  inputSignal.push_back(Eigen::VectorXd::Zero(sizes[0]));
  //listDetsToTrain.push_back(fullBasis.getDetByIndex(0));
  for (int layer=0; layer < numLayersBiasesWeights; ++layer){
    activations.push_back(Eigen::VectorXd::Zero(sizes[layer+1]));
    inputSignal.push_back(Eigen::VectorXd::Zero(sizes[layer+1]));
    
    gFactorBiases.push_back(Eigen::VectorXd::Ones(sizes[layer+1]));
    gFactorBiasesPrev.push_back(Eigen::VectorXd::Ones(sizes[layer+1]));

    gFactorWeights.push_back(Eigen::MatrixXd::Ones(sizes[layer+1], sizes[layer]));
    gFactorWeightsPrev.push_back(Eigen::MatrixXd::Ones(sizes[layer+1], sizes[layer]));
  }
  //nablaWeightsPrev = nablaWeights;
  //nablaBiasesPrev = nablaBiases;

  // Assign the output state and the evaluator to 0
  outputState = State();
}

void NeuralNetwork::train(std::vector<detType> const &listDetsToTrain, double eta, double epsilon){
// The coefficients are stored in scope of the train method and then stored into the state
  //nablaWeightsPrev = nablaWeights;
  //nablaBiasesPrev = nablaBiases;
  int numDets(listDetsToTrain.size());
  std::vector<std::vector<Eigen::VectorXd>> inputSignalEpochs;
  std::vector<std::vector<Eigen::VectorXd>> activationsEpochs;
  std::vector<detType> listDetsToTrainFiltered;
  std::vector<coeffType > outputCs;
  std::vector<std::vector<coeffType>> coupledOutputCsEpochs;
  std::vector<coeffType> coupledOutputCs;
  std::vector<std::vector<detType>> coupledDetsEpochs;
  std::vector<detType> coupledDets;
  int numLayersNeuron(sizes.size());
  int numLayersBiasesWeights(sizes.size()-1);
  //reset nablaWeights and nablaBiases to 0 values;
  nablaPara *= 0.;
  //for (int layer=0; layer < numLayersBiasesWeights; ++layer){
  //  nablaBiases[layer] *= 0.; 
  //  nablaWeights[layer] *=0.;
  //}

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
}

void NeuralNetwork::updateParameters(int method){
  //method corresponds to
  // 0: Stochastic gradiend desend
  // 1: Stochastic reconfiguration
  Eigen::VectorXd generlisedForce=backPropagate(inputSignalEpochs, activationsEpochs);
  if (method == 0){
    NNP = sl.update(NNP,generlisedForce);
    //update weights and biases
    //0th layer is the input layer,
    //the biases should not change.
    //weights have only 0 -- numLayers-2 layers.
 
    //stochastic gradient descend 
    /*
    for (int layer=0; layer < numLayersBiasesWeights; ++layer){
      if(momentum){
        gFactorBiases[layer] = ((nablaBiases[layer].array() 
                              * nablaBiasesPrev[layer].array()) > 1e-8).select(
                              (gFactorBiasesPrev[layer].array()+0.05).matrix(), gFactorBiasesPrev[layer]*0.95);
        gFactorWeights[layer] = ((nablaWeights[layer].array()
                         * nablaWeightsPrev[layer].array()) > 1e-8).select(
                           (gFactorWeightsPrev[layer].array()+0.05).matrix(), gFactorWeightsPrev[layer]* 0.95);
        nablaBiases[layer] = (nablaBiases[layer].array() * gFactorBiases[layer].array()).matrix();
        nablaWeights[layer] = (nablaWeights[layer].array() * gFactorWeights[layer].array()).matrix();
      }
      //nablaBiases[layer] = -eta/numDets  * nablaBiases[layer];
      //nablaWeights[layer] = -eta/numDets  * nablaWeights[layer]; 
      nablaBiases[layer] = -eta  * nablaBiases[layer];
      nablaWeights[layer] = -eta  * nablaWeights[layer]; 
      if (momentum){
        nablaBiases[layer] += momentumDamping * nablaBiasesPrev[layer];
        nablaWeights[layer] += momentumDamping * nablaWeightsPrev[layer]; 
      }
      biases[layer] += nablaBiases[layer];
      weights[layer] += nablaWeights[layer];
    }
  */
  }
  //Stochastic reconfiguration
  else {
    Eigen::MatrixXd dCdw=backPropagateSR(inputSignalEpochs, activationsEpochs);
    Eigen::VectorXd ci = Eigen::Map<Eigen::VectorXd>(&(outputState.getAllCoeff()[0]),outputState.size());
    NNP = sl.update(NNP,generlisedForce,ci,dCdw);
     
    //store dCdw
    //
    
    //construct matrix s_{k.k'}.
    //starting with vector <O_k>. Need outputState (global class member variable)

  }
  

  //double probAmp(0.);
  //double max(0);
  //for (size_t i=0; i < outputCs.size(); ++i){
  //  probAmp = std::norm(outputCs[i]); //probility amplitude
  //  if (fabs(probAmp) -fabs(max) > 1e-8) max = probAmp;
  //}
  //std::random_device rngd;
  //double const normalizerd = static_cast<double>(rngd.max());
  //std::vector<detType> seeds;
  //detType seedPrev;
  //seeds.push_back(listDetsToTrain[0]);
  //double cPrev = std::norm(outputCs[0]);
  //for (size_t i=0; i < outputCs.size(); ++i){
  //  probAmp = std::norm(outputCs[i]); //probility amplitude
  //  double prandom=rngd()/normalizerd;
  //  if (fabs(probAmp)-epsilon*fabs(max) > 1e-8){
  //  //if (prandom-probAmp<-1e-8){
  //    seeds.push_back(listDetsToTrainFiltered[i]);
  //    cPrev = std::norm(outputCs[i]);
  //  }
  //}
  //return seeds;
}

Eigen::VectorXd NeuralNetwork::feedForward(detType const& det) const{
  int numStates=det.size();
  int numLayersNeuron(sizes.size());
  for (int state=0; state<numStates; ++state){
    //std::cout << "det[" << state << "]=" << det[state] << std::endl;
    activations[0][state] = det[state]?1.0:0.0;
    inputSignal[0][state] = det[state]?1.0:0.0;
  }

  for (int layer=1; layer < numLayersNeuron; ++layer){
    //Here the 0th layer of weights correspond to the connections between
    //the 0th and 1st layer. Here layer refers to the Neuron layer. When use 
    //it to refer the Biases and weights' layer, we need conversion.
    //map the para into layered structure.
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
       std::vector<std::vector<Eigen::VectorXd>> const &activationsEpochs,
     ){
  int numDets = outputState.size();
  int numLayersNeuron(sizes.size());
  int numLayersBiasesWeights(sizes.size()-1);
  // catch here
  //everytime the backPropagate is called, we should reset nabla* to zero.
  nablaNNP *= 0.;
  //for (int layer=0; layer < numLayersBiasesWeights; ++layer){
  //  nablaBiases[layer] *= 0.; 
  //  nablaWeights[layer] *=0.;
  //}
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

      //check if it right
      deltaThisLayer = deltaThisLayer.array()
         * inputSignalEpochs[epoch][layer+1].unaryExpr(&Tanh_prime).array();
      //deltaThisLayer = deltaAsDiagonal * 
      // VectorXd::Ones(deltaThisLayer.size());
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

Eigen::MatrixXd NeuralNetwork::backPropagateSR(
       std::vector<std::vector<Eigen::VectorXd>> const &inputSignalEpochs,
       std::vector<std::vector<Eigen::VectorXd>> const &activationsEpochs,
     ){
  int numDets = outputState.size();
  int numLayersNeuron(sizes.size());
  int numLayersBiasesWeights(sizes.size()-1);
  // catch here
  //everytime the backPropagate is called, we should reset nabla* to zero.
  //for (int layer=0; layer < numLayersBiasesWeights; ++layer){
  //  nablaBiases[layer] *= 0.; 
  //  nablaWeights[layer] *=0.;
  //}
  std::vector<Eigen::VectorXd> dEdC(numDets,Eigen::VectorXd::ones(2));
  //create w_k|--C_i matrix
  Eigen::MatrixXd dCdw(numPara,numDets);
  for (int epoch(0); epoch < numDets; ++epoch){
    //reset nablaNNP to 0
    nablaNNP *= 0.;
    //for (int layer=0; layer < numLayersBiasesWeights; ++layer){
    //  nablaBiases[layer] *= 0.; 
    //  nablaWeights[layer] *=0.;
    //}
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

      //check if it right
      deltaThisLayer = deltaThisLayer.array()
         * inputSignalEpochs[epoch][layer+1].unaryExpr(&Tanh_prime).array();
      //deltaThisLayer = deltaAsDiagonal * 
      // VectorXd::Ones(deltaThisLayer.size());
      deltaAllLayers.push_back(deltaThisLayer);
      //combine all nablaWeights and nablaBiases to a single vector
      nablaBiases[layer] = deltaThisLayer;
      //get a weight matrix. \partial C/\partial w^{l}_{jk} = a^{l-1}_k \delta_j^l 
      //the layer here refers to the lth layer of Biases and weights, so for
      //activation layer refers to the l-1th layer.
      //reshape matrix into a vector
      nablaWeights[layer] = deltaThisLayer * activationsEpochs[epoch][layer].transpose();
      //reshape nablaWeights into a vector
    }
    //fill up the dCdw matrix
    dCdw.col(epoch) << nablaNNP; 
  }
  return dCdw;
}

void preTrain(NeuralNetwork &network, State const &target, double trainRate, double epsilon){
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
		network.train(list, trainRate, epsilon);
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
