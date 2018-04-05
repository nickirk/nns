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
  int sizeFilter;
  //numNrn here refers to length of the filter, so in total
  //there should be numNrn*depthLastLayer parameters per filter.
  int numNrn;
  int stride;
  std::vector<Eigen::VectorXd> convolve(Eigen::VectorXd const& input);
  //z=w*input+bias, which is an intermediate variable but is needed
  //during the backpropagation step
  std::vector<Eigen::VectorXd> z;
  std::vector<Eigen::Map<Eigen::VectorXd>> biases;
  std::vector<Eigen::Map<Eigen::MatrixXd>> weights;

}


#endif
