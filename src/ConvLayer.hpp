/*
 * ConvLayer.hpp
 * Created on 27.3.2018
 * Author: Ke Liao 
 */
#ifndef ConvLayer_DEFINED
#define ConvLayer_DEFINED

class ConvLayer: public Layer {
public:
  ConvLayer(
            std::vector<Eigen::VectorXd> const &inputs_, 
            double &(actFunc_)(double), int numFilters_, 
            int size_, int stride_
            ) 
  virtual ~ConvLayer();
  void processSignal();
  void mapPara(double *adNNP, int &startPoint);
  void backProp();
  int numPara;
private:
  int numFilters;
  int depthFilter;
  int lengthFilter;
  // numNrn here refers to length of the filter, so in total
  // there should be numNrn*depthLastLayer parameters per filter.
  int numNrn;
  int stride;
  std::vector<Eigen::VectorXd> convolve();
  // z=w*input+bias, which is an intermediate variable but is needed
  // during the backpropagation step
  std::vector<Eigen::VectorXd> z;
  // each filter (consists of several chanels) shares the same 
  // bias. So the dimension will be L vector of double
  std::vector<Eigen::Map<Eigen::VectorXd>> biases;
  // the weights have one-to-one correspondance to neurons. 
  // in convNet, each filter has the same depth as that of the 
  // previous output signal. So we define it as a Lxlxsize
  std::vector<std::vector<Eigen::Map<Eigen::MatrixXd>>> weights;

}


#endif
