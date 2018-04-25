/*
 * SpinConfig.hpp
 *
 *  Created on: Apr 24, 2018
 *      Author: guther
 */

#ifndef SRC_UTILITIES_SPINCONFIG_HPP_
#define SRC_UTILITIES_SPINCONFIG_HPP_


// Storing how many electrons have which ms
class SpinConfig {
public:
	// for now, only differ between ms>0 and ms<0
	SpinConfig(int up_, int down_, int numStates_):up(up_),down(down_),numStates(numStates_){};
	virtual ~SpinConfig(){};
	// give an integer ms and get the number of electrons with that spin
	int operator()(int ms) const{
		if(ms == 1) return up;
		return down;
	}
	// total number of spin orbs
	int numSpinOrbs() const{return numStates;}
private:
	// internal stuff
	int up, down, numStates;
};

#endif /* SRC_UTILITIES_SPINCONFIG_HPP_ */
