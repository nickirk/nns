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
// get the sign of the jordan-wigner string from annihilatorIndex to creatorIndex applied to a
int JWStringSign(detType const &a, int annihilatorIndex, int creatorIndex);
// Get the occupied orbitals
std::vector<int> getOccupiedPositions(detType const &det);
// Read the vector of bools as an integer (i.e. reinterpret cast-style)
int verbatimCast(detType const & det);

// for testing purposes: get the number of orbitals in which a,b differ
int getExcitLvl(detType const &a, detType const &b);
// useful for matrix element generation: get the orbs in which a and b differ
void getExcitation(detType const &a, detType const &b, std::vector<int> &excitations,
		std::vector<int> &holes, std::vector<int> &same);

// binary operators for determinants
bool operator == (detType const &a, detType const &b);
bool operator < (detType const &lhs, detType const &rhs);
bool compare_det(detType const &lhs, detType const &rhs);

// output function for debugging
void printDet(detType const &out);

}
#endif
