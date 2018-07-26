/*
 * errors.hpp
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */

#ifndef ERRORS_HPP_
#define ERRORS_HPP_
#include <string>
// classes for error handling. This is the way to go for
// exception handling

namespace networkVMC{

class OutOfRangeError{
public:
	OutOfRangeError(int pos_):pos(pos_){};
	int pos;
};

class InvalidAnnihilation{
public:
	InvalidAnnihilation(int pos_):pos(pos_){};
	int pos;
};

class InvalidCreation{
public:
	InvalidCreation(int pos_):pos(pos_){};
	int pos;
};

class SizeMismatchError{
public:
	SizeMismatchError(int a,int b):sizeA(a),sizeB(b){};
	int sizeA, sizeB;
};

class ActFuncDoNotExist{
  public:
    ActFuncDoNotExist(std::string actFunc_):actFunc(actFunc_){};
    std::string actFunc;
};

class InvalidParameterPassed{
public:
	InvalidParameterPassed(){};
};

class EmptyNetworkError{
public:
	EmptyNetworkError(){};
};

class NoExcitationFound{
public:
	NoExcitationFound(int lvl_):lvl(lvl_){};
	int lvl;
};

class UnconvergedEigenproblem{
public:
	UnconvergedEigenproblem(){};
};

class InvalidDeterminantError{
public:
	InvalidDeterminantError(){};
};

}

#endif /* ERRORS_HPP_ */
