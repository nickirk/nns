/*
 * Nns.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <limits>
#include <math.h>
#include <string>
#include <iomanip>
#include "Nns.hpp"


int main(){
  Hamiltonian ham(10, 0.95, 0.2, 0.4);
  std::cout << ham.entry << std::endl;
}
