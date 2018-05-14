#include <random>
#include <iostream>

int main(){
  std::random_device rd;     // only used once to initialise (seed) engine
  std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
  std::uniform_real_distribution<double> uni(0,1); // guaranteed unbiased
  for (int i(0); i<10; ++i){
    std::cout << uni(rd) << std::endl;
  }
}
