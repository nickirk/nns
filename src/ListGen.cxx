/*
 * ListGen.cxx
 *
 *  Created on: Jan 23, 2018
 *      Author: guther
 */

#include "Sampler.hpp"
#include "ListGen.hpp"

ListGen::ListGen(Hamiltonian const &H_, Basis const &fullBasis_, int numDets_, detType const &HF, NeuralNetwork const &NNW_):
	Sampler(H_,fullBasis_,numDets_,HF),NNW(NNW_),pos(0){
	diffuseList=std::vector<detType >(numDets_,HF);
}

ListGen::~ListGen() {
}

void ListGen::iterate(coeffType &cI, detType &dI) const{
	// Fetch the next entry from the pre-arranged list
	dI = getDet();
	// Get its coefficient
	cI = NNW.getCoeff(dI);
	// Dont forget to cache the state
	NNW.cacheNetworkState();
}

detType ListGen::getDet() const{
	pos += 1;
	// cycle through the list, if we reach the end, start from the beginning
	if(pos > diffuseList.size()) pos = 1;
	return diffuseList[pos-1];
}

void ListGen::diffuse(std::vector<detType> &list, std::vector<int> const& spinConfig){
 list.clear();
 detType buf;
 //std::cout << "before adding random dets=" << list.size() << std::endl;
 while (list.size() < numDets){
   buf = getRandomDeterminant(spinConfig);
   list.push_back(buf);
 }
 //std::cout << "after adding random dets=" << list.size() << std::endl;
 coeffType c_i;
 coeffType c_j;
 Eigen::VectorXd lastLayerActivation;
 std::random_device rngd;
 double const normalizerd = static_cast<double>(rngd.max());
 double prandom=rngd()/normalizerd;
 double probUnbias(0.);
 for (size_t i(0); i<list.size(); ++i){
   prandom=rngd()/normalizerd;
   lastLayerActivation=NNW.feedForward(list[i]);
   c_i=coeffType(lastLayerActivation[0], lastLayerActivation[1]);
   buf = getRandomCoupledState(list[i],probUnbias);
   //buf = getRandomDeterminant(spinConfig);
   lastLayerActivation = NNW.feedForward(buf);
   c_j=coeffType(lastLayerActivation[0], lastLayerActivation[1]);
   //getRandomCoupledState(buf,probUnbias);
   if (prandom - std::pow(std::norm(c_j),2)/std::pow(std::norm(c_i),2)<-1e-8){
     list[i]=buf;
   }
 }
 //std::cout << "after diffuse=" << list.size() << std::endl;
 removeDuplicate(list);
 diffuseList = list;
 pos = 0;
 //std::cout << "after remove dups=" << list.size() << std::endl;
}
/*
void Sampler::generateList(std::vector<detType > &list) const{
  // this just repeatedly gets new random states and adds them to the list
  // if the number of target states is not an integer multiple of the number of
  // reference states, we have to round up the number of target states
  //int numComp = numDets;
  //list = std::vector<detType >(numDets);
  list.clear();
  detType buf;
  detType bufPrev;
  int numRef=cDet.size();
  int random_integer(0);
  //if(numDets%cDet.size()!=0) numComp = (numDets/cDet.size()+1)*cDet.size();
  //std::cout << "size of numComp= " << numComp << std::endl;
  std::random_device rd;     // only used once to initialise (seed) engine
  std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
  std::uniform_int_distribution<int> uni(0,numRef-1); // guaranteed unbiased
  std::random_device rngd;
  double const normalizerd = static_cast<double>(rngd.max());
  if (numRef < numDets) {
    for (int i(0); i<numRef; ++i){
      list.push_back(cDet[i]);
    }
    //std::random_device rng;     // only used once to initialise (seed) engine
    //double const normalization=static_cast<double>(rng.max());
    //random_integer=static_cast<int>(rng()/normalization*(numRef-1));
    //for (int i(numRef); i < numDets; ++i){
   double prandom=rngd()/normalizerd;
   if (prandom-0.8 < -1e-8){
    random_integer = uni(rng);
    buf = cDet[random_integer];
   }
   else{
    buf = getRandomDeterminant(3,3,16);
   }
    bufPrev = buf;
    while (list.size() < numDets){
      prandom=rngd()/normalizerd;
      if (prandom-1. < -1e-8){
      //if (true){
        //for (int depth(0);  depth < 4; depth++){
          while (verbatimCast(bufPrev)==verbatimCast(buf) ){
            buf = bufPrev;
            buf = getRandomConnection(buf);
          }
          list.push_back(buf);
          bufPrev = buf;
        //}
      }
      else {
        buf = getRandomDeterminant(3,3,16);
        list.push_back(buf);
      }
    }
  }
  else if (numRef == numDets) list = cDet;
}
*/

