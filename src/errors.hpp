/*
 * errors.hpp
 *
 *  Created on: Aug 14, 2017
 *      Author: guther
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

#endif /* ERRORS_HPP_ */
