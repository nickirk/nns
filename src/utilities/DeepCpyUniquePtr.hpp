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

/**
 * \class DeepCpyUniquePtr
 * \brief A smart pointer of clonable objects
 * \tparam T type of the managed object
 *
 * This is a unique pointer for polymorphic objects featuring a virtual constructor
 * (clone member function). It can be copy-constructed and assigned, in this case, a deep copy of
 * the managed resource is created.
 */
template <typename T>
class DeepCpyUniquePtr {
  public:
	/// Default initialization: resource is set to nullptr
	DeepCpyUniquePtr():resource(nullptr){};
	/**
	 * \brief Resource initialization: take ownership of an object
	 * \param source pointer to the object to manage. Ownership is taken
	 */
	DeepCpyUniquePtr(T *source):resource(source){};
	// also allow for construction from unique_ptr, to circumvent usage of raw pointers

	/**
	 * \brief Resource initialization: take ownership of an object
	 * \param source transfers ownership from the source to *this
	 */
	DeepCpyUniquePtr(std::unique_ptr<T> &source):resource(source.release()){};

// Copy/Move operations (constructors/assignment operators)

//---------------------------------------------------------------------------//

	/**
	 * \brief deep copy-ctor
	 * \param[in] source DeepCpyUniquePtr to be copied
	 */
	DeepCpyUniquePtr(DeepCpyUniquePtr<T> const &source):resource(nullptr){
		if(source.resource != nullptr) resource = source.resource->clone();
	}

//---------------------------------------------------------------------------//

	/**
	 * \brief move-ctor (works intuitively)
	 * \param source DeepCpyUniquePtr from which the managed resource shall be moved
	 */
	DeepCpyUniquePtr(DeepCpyUniquePtr<T> &&source):resource(nullptr){
		swap(*this,source);
	}

//---------------------------------------------------------------------------//

	/**
	 * \brief deep copy/move assignment
	 * \param[in] source DeepCpyUniquePtr to be copied/moved
	 */
	DeepCpyUniquePtr& operator=(DeepCpyUniquePtr<T> source){
		// we create a deep copy, in particular, T has to be clonable
	    swap(*this,source);
		return *this;
	}

//---------------------------------------------------------------------------//

	// Also, we need the pointer functionality
	/**
	 * \brief const operator to de-reference
	 * \return const reference to the managed object
	 */
	T const& operator*()const {return *resource;}
	/**
	 * \brief const operator to call member functions
	 * \return const pointer to the managed object and recursively calls operator->()
	 */
	T const* operator->()const {return resource;}

	/**
	 * \brief operator to de-reference
	 * \return reference to the managed object
	 */
	T& operator*() {return const_cast<T&>(
			*static_cast<DeepCpyUniquePtr<T> const&> (*this));
	}
	/**
	 * \brief const operator to call member functions
	 * \return pointer to the managed object and recursively calls operator->()
	 */
	T* operator->() {return resource;}

	// RAII
	/// deletes the managed resource
	virtual ~DeepCpyUniquePtr(){delete resource;}


	//---------------------------------------------------------------------------//

	/// \brief swaps the content of a and b
	friend void swap(DeepCpyUniquePtr<T> &a, DeepCpyUniquePtr<T> &b){
		std::swap(a.resource, b.resource);
	}
private:
	/// \brief managed resource
	T* resource;
};

} /* namespace networkVMC */

#endif /* SRC_UTILITIES_DEEPCPYUNIQUEPTR_HPP_ */
