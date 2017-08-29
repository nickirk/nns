#include <random>
#include <iostream>

int main(){
  std::random_device rd;     // only used once to initialise (seed) engine
  std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
  std::uniform_int_distribution<int> uni(0,10); // guaranteed unbiased
  int random_integer;
  for (int i(0); i<10; ++i){
    random_integer = uni(rng);
    std::cout << random_integer << std::endl;
  }
}
