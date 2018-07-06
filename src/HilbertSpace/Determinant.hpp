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

#include "../utilities/TypeDefine.hpp"
#include "../utilities/Errors.hpp"

namespace networkVMC{

// these should not be members of some class, better use generic functions
void create(detType &det, int pos);
// Remove an electron at orbital pos
void annihilate(detType &det, int pos);
// Get the occupied orbitals
std::vector<int> getOccupiedPositions(detType const &det);
// Read the vector of bools as an integer (i.e. reinterpret cast-style)
int verbatimCast(detType const & det);

// for testing purposes: get the number of orbitals in which a,b differ
int getExcitLvl(detType const &a, detType const &b);

// binary operators for determinants
bool operator == (detType const &a, detType const &b);
bool operator < (detType const &lhs, detType const &rhs);
bool compare_det(detType const &lhs, detType const &rhs);

// output function for debugging
void printDet(detType const &out);

}
#endif
