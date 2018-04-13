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
#include "utilities/Errors.hpp"
#include "TypeDefine.hpp"
// these should not be members of some class, better use generic functions
void create(detType &det, int pos);
void annihilate(detType &det, int pos);
std::vector<int> getOccupiedPositions(detType const &det);
int verbatimCast(detType const & det);
bool operator == (detType const &a, detType const &b);
bool operator < (detType const &lhs, detType const &rhs);
bool compare_det(detType const &lhs, detType const &rhs);
#endif
