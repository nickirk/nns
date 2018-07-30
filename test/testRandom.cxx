#include "../src/utilities/RNGWrapper.hpp"
#include <iostream>

int main(){
  networkVMC::RNGWrapper rng;
  for (int i(0); i<10; ++i){
    std::cout << rng() << std::endl;
  }
}
