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

namespace networkVMC{

// these should not be members of some class, better use generic functions
void create(detType &det, int pos); // can throw InvalidCreation
// Remove an electron at orbital pos
void annihilate(detType &det, int pos); // can throw InvalidAnnihilation
// excite from orbtial i to j
detType excite(detType const &source, int i, int j); // can throw InvalidCreation/Annihilation and OutOfRangeError
// get the exponent of the jordan-wigner string sign from annihilatorIndex to creatorIndex applied to a
int JWStringLength(detType const &a, int annihilatorIndex, int creatorIndex);
// sign induces by applying a_crt^dagger a_ann on a
inline int excitationSign(detType const &a, int ann, int crt){
	return -2*(JWStringLength(a,ann,crt)%2)+1;
}
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

// remove duplicates from a list of determinants
void removeDuplicate(std::vector<detType> &list);

}
#endif
