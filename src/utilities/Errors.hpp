/*
 * errors.hpp
 * Copyright (c) 2017
 * Author: Ke Liao and Kai Guther
 * All rights reserved
 */

#ifndef ERRORS_HPP_
#define ERRORS_HPP_

// classes for error handling. This is the way to go for
// exception handling

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

//class ActFuncDoNotExist{
//  public:
//    ActFuncDoNotExist(std::string actFunc_):actFunc(actFunc_);
//    std::string actFunc;
//};

#endif /* ERRORS_HPP_ */
