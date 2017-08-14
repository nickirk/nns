#ifndef Determinant_DEFINED
#define Determinant_DEFINED

#include <stdio.h>
#include <vector>
class Determinant{
   public:
     Determinant();
     Determinant(int size_);
     Determinant(Determinant const &determinant_);
     void operator = (Determinant const &determinant_);
     void annihilate(int pos);
     void create(int pos); 
     std::vector<int> getOccupiedPositions() const;
     int getSize() const;
     int intCast() const;
     bool operator <(Determinant const &det_) ;
   private:
     int size;
     std::vector<bool> det;
};


#endif
