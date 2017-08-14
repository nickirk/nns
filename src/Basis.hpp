#ifndef Basis_DEFINED
#define Basis_DEFINED

#include "Determinant.hpp"
#include "detType.hpp"
//
// type of the determiants


//class BasisIndexRef{
//  public:
//    BasisIndexRef(): index(0){}
//    BasisIndexRef(Determinant det_, int index_): det(det_), index(index_){}
//    BasisIndexRef(BasisIndexRef const & basisIndexRef): det(basisIndexRef.det),
//      index(basisIndexRef.index){}
//    static bool sortByIndex(BasisIndexRef const& basisindexref1,
//		  BasisIndexRef const& basisindexref2){
//      return basisindexref1.index < basisindexref2.index;
//    };
//    Determinant det;
//    int index;
//};

class Basis{
  public:
    Basis(int numEle_, int numOrb_);
    int getSize();
    detType getDetByIndex(int index);
    int getIndexByDet(detType const & det_);
  private:
    int numEle;
    int numOrb;
    int size;
    int indexOfDet;
    std::vector<int> listOfOrbNum;
    std::vector<int> combination;
    int calcSize(int numOrb_, int numEle_);
    void createBasisDet(int offset, int numEle_);
    // determinants need to be accessed often, no need to add some
    // overhead by using an extra class here, better use an alias
    std::vector<detType > basis;
    std::vector<int> indexBasis;
};

#endif
