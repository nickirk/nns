/*
 * testStateSort.cxx
 *
 *  Created on: May 4, 2018
 *      Author: guther
 */

#include <random>
#include <iostream>
#include "../src/utilities/State.hpp"

using namespace networkVMC;

int main(){
  State psi;
  psi.resize(10);
  std::random_device rnd;
  double normalizer = rnd.max();
  std::cout<<"Unsorted state coefficients"<<std::endl;
  for(std::size_t i = 0; i<psi.size(); ++i){
	  psi.coeff(i) = rnd()/normalizer;
	  std::cout<<psi.coeff(i)<<std::endl;
  }
  psi.sortCoeffs();
  std::cout<<"Sorted state coefficients"<<std::endl;
  for(std::size_t i = 0; i<psi.size(); ++i){
	  std::cout<<psi.coeff(i)<<std::endl;
  }
}

