#ifndef Determinant_DEFINED
#define Determinant_DEFINED

#include <stdio.h>
#include <vector>
class Determinant{
   public:
     Determinant(int size_);
     Determinant(Determinant &determinant_);
     void annihilate(int pos);
     void create(int pos); 
     std::vector<int> getOccupiedPositions() const;
     int getSize() const;
     int intCast() const;
   private:
     int size;
     std::vector<bool> det;
};


#endif
