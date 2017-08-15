/*
 * Hamiltonian.cxx
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */

#ifndef Determinant_DEFINED
#define Determinant_DEFINED

#include <stdio.h>
#include <vector>
#include "detType.hpp"

// deprecated
class Determinant{
   public:
     Determinant();
     Determinant(int size_);
     Determinant(Determinant const &determinant_);
     void operator = (Determinant const &determinant_);
     void annihilate(int pos);
     void create(int pos); 
     std::vector<int> getOccupiedPositions() const;
     int getSize() const;
     int intCast() const;
     bool operator == (Determinant const &det_) const;
   private:
     int size;
     detType det;
};


// these should not be members of some class, better use generic functions
void create(detType &det, int pos);
void annihilate(detType &det, int pos);
// this operator likely never occurs, better use a lambda when needed
bool operator<(detType const &a, detType const &b);

#endif
