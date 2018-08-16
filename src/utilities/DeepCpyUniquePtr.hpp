/*
 * DeepCpyUniquePtr.hpp
 *
 *  Created on: Jun 20, 2018
 *      Author: guther
 */

#ifndef SRC_UTILITIES_DEEPCPYUNIQUEPTR_HPP_
#define SRC_UTILITIES_DEEPCPYUNIQUEPTR_HPP_

#include <memory>

namespace networkVMC {

// a smart pointer of clonable objects
template<typename T>
class DeepCpyUniquePtr {
public:
	DeepCpyUniquePtr():resource(nullptr){};
	DeepCpyUniquePtr(T *source):resource(source){};
	// also allow for construction from unique_ptr, to circumvent usage of raw pointers
	// this transfers ownership from the unique_ptr to *this
	DeepCpyUniquePtr(std::unique_ptr<T> &source):resource(source.release()){};

// Copy/Move operations (constructors/assignment operators)

//---------------------------------------------------------------------------//

	DeepCpyUniquePtr(DeepCpyUniquePtr const &source):resource(nullptr){
		if(source.resource != nullptr){
			makeClone(source);
		}
	}

//---------------------------------------------------------------------------//

	DeepCpyUniquePtr(DeepCpyUniquePtr &&source):resource(nullptr){
		std::swap(source.resource,resource);
	}

//---------------------------------------------------------------------------//

	DeepCpyUniquePtr& operator=(DeepCpyUniquePtr const &source){
		// we create a deep copy, in particular, T has to be clonable
		makeClone(source);
		return *this;
	}

//---------------------------------------------------------------------------//

	DeepCpyUniquePtr& operator=(DeepCpyUniquePtr &&source){
		// swap with source - deallocation will be done automatically in
		// destruction of source
		std::swap(source.resource,resource);
		return *this;
	}

//---------------------------------------------------------------------------//

	// Also, we need the pointer functionality
	T const& operator*()const {return *resource;}
	T const* operator->()const {return resource;}

	T& operator*() {return const_cast<T&>(
			*static_cast<DeepCpyUniquePtr<T> const&> (*this));
	}

	T* operator->() {return resource;}

	// RAII
	virtual ~DeepCpyUniquePtr(){delete resource;}
private:
	void makeClone(DeepCpyUniquePtr const& source){
		// make a deep copy
		resource = source.resource->clone();
	}

	T* resource;
};

} /* namespace networkVMC */

#endif /* SRC_UTILITIES_DEEPCPYUNIQUEPTR_HPP_ */
