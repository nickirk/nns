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

class outOfRangeError{
public:
	outOfRangeError(int pos_):pos(pos_){};
	int pos;
};

class invalidAnnihilation{
public:
	invalidAnnihilation(int pos_):pos(pos_){};
	int pos;
};

class invalidCreation{
public:
	invalidCreation(int pos_):pos(pos_){};
	int pos;
};

class sizeMismatchError{
public:
	sizeMismatchError(int a,int b):sizeA(a),sizeB(b){};
	int sizeA, sizeB;
};

#endif /* ERRORS_HPP_ */
