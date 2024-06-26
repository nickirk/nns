/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARDFMat.h
   Matrix template that generates a dense matrix from a file.

   ARPACK authors:
      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/


#ifndef ARDFMAT_H
#define ARDFMAT_H

#include <cstddef>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <string>

#include "../../../lib/arpackpp/include/arch.h"
#include "../../../lib/arpackpp/include/arerror.h"


template<class ARTYPE>
class ARdfMatrix {

 private:

  // const int linelength = 256;

  std::string   datafile;  // Filename.
  std::ifstream file;      // File handler.
  int      m;         // Number of rows.
  int      n;         // Number of columns.
  int      blksize;   // Size of each matrix block that is read at once.
  int      block;     // Index of the matrix block that is to be read.
  int      nblocks;   // Number of blocks the matrix contain. 
  int      first;     // First row/column stored in val.
  int      strows;    // Number of rows actually stored in val.
  int      stcols;    // Number of columns actually stored in val.
  int      headsize;  // Number of lines in the heading part of the file
                      // (including the line that contains the matrix size).
  bool     roword;    // A variable that indicates if the data will be read
                      // using a row-major or a column-major ordering.
  ARTYPE*  val;       // Numerical values of matrix entries.
			
  void ConvertDouble(char* num);

  void SetComplexPointers(char* num, char* &realp, char* &imagp);
			  			  
  bool ReadEntry(std::ifstream& file, double& val);

  bool ReadEntry(std::ifstream& file, float& val);

  bool ReadEntry(std::ifstream& file, arcomplex<double>& val);

  bool ReadEntry(std::ifstream& file, arcomplex<float>& val);
  
 public:
  
  bool IsDefined() const { return (m!=0); }
  
  bool IsOutOfCore() const { 
    return ((m!=0) && ((roword && (blksize<m)) || (blksize<n))); 
  }

  bool IsRowOrdered() const { return roword; }

  std::string Filename() const { return datafile; }

  void Rewind();

  int BlockSize() const { return blksize; }

  int FirstIndex() const { return first; }

  int NBlocks() const { return nblocks; }

  int ColsInMemory() const { return stcols; }
  
  int RowsInMemory() const { return strows; }
  
  int NRows() const { return m; }

  int NCols() const { return n; }
  
  ARTYPE* Entries() const { return val; }

  void ReadBlock();
  // Function that reads a block of blksize rows/columns of the matrix.

  void Define(const std::string& filename, int blksizep = 0);
  // Function that reads the matrix dimension. Define also read all
  // of the matrix elements when blocksize = 0.
  
  ARdfMatrix();
  // Short constructor.

  ARdfMatrix(const std::string& filename, int blksizep = 0);
  // Long constructor.

  ~ARdfMatrix();
  // Destructor.

}; // Class ARdfMatrix.


// ------------------------------------------------------------------------ //
// ARdfMatrix member functions definition.                                  //
// ------------------------------------------------------------------------ //


template<class ARTYPE>
inline void ARdfMatrix<ARTYPE>::ConvertDouble(char* num) 
{

  char* pd;

  pd = strchr((char*)num,'D');
  if (pd) *pd = 'E';
  pd = strchr((char*)num,'d');
  if (pd) *pd = 'E';

} // ConvertDouble.


template<class ARTYPE>
inline void ARdfMatrix<ARTYPE>::
SetComplexPointers(char* num, char* &realp, char* &imagp)
{

  realp = num;
  while (*realp == ' ') realp++;
  imagp = realp;
  while (*imagp != ' ') imagp++;  

} // SetComplexPointers.


template<class ARTYPE>
inline bool ARdfMatrix<ARTYPE>::ReadEntry(std::ifstream& file, double& val)
{

  char num[LINELEN];
  char c;

  if (file.get((char*)num,LINELEN,'\n')) {
    file.get(c);
    ConvertDouble((char*)num);
    val = atof((char*)num);
    return true;
  }
  else {
    return false;
  }

} // ReadEntry (double).


template<class ARTYPE>
inline bool ARdfMatrix<ARTYPE>::ReadEntry(std::ifstream& file, float& val)
{

  double dval;
  bool   ret;
  
  ret = ReadEntry(file, dval);
  val = (float)dval;
  return ret;

} // ReadEntry (float).


template<class ARTYPE>
inline bool ARdfMatrix<ARTYPE>::
ReadEntry(std::ifstream& file, arcomplex<double>& val)
{

  char  num[LINELEN];
  char  c;
  char  *realp, *imagp;

  if (file.get((char*)num,LINELEN,'\n')) { 
    file.get(c);
    SetComplexPointers((char*)num, realp, imagp);
    ConvertDouble((char*)realp);
    ConvertDouble((char*)imagp);
    val = arcomplex<double>(atof((char*)realp), atof((char*)imagp));
    return true;
  }
  else {
    return false;
  }

} // ReadEntry (arcomplex<double>).


template<class ARTYPE>
inline bool ARdfMatrix<ARTYPE>::
ReadEntry(std::ifstream& file, arcomplex<float>& val)
{

  char  num[LINELEN];
  char  c;
  char  *realp, *imagp;

  if (file.get((char*)num,LINELEN,'\n')) { 
    file.get(c);
    SetComplexPointers((char*)num, realp, imagp);
    ConvertDouble((char*)realp);
    ConvertDouble((char*)imagp);
    val = arcomplex<float>(atof((char*)realp), atof((char*)imagp));
    return true;
  }
  else {
    return false;
  }

} // ReadEntry (arcomplex<float>).


template<class ARTYPE>
void ARdfMatrix<ARTYPE>::Rewind() 
{ 

  char data[LINELEN];
  char c;

  file.seekg(0);
  block  = 0; 
  first  = 0;
  strows = 0;
  stcols = 0;

  // Skipping the header. 

  for (int i=0; i<headsize; i++) {
    file.get((char*)data,LINELEN,'\n');
    file.get(c);
  }

} // Rewind.


template<class ARTYPE>
void ARdfMatrix<ARTYPE>::ReadBlock()
{

  int    i, j, last;
  ARTYPE value;

  // Repositioning the file pointer if block == 0.

  if (block == 0) Rewind();

  // Reading a block.

  first = (block++)*blksize; // First row/column to be read.
  last  = first+blksize;     // First row/column of the next block.

  if (roword) {

    // Adjusting last if we are going to read the last block.

    if (last > m) {
      last  = m;
      block = 0;
    }
    last  -= first;
    strows = last;
    stcols = n;

    // Reading matrix data.

    for (i=0; i<last; i++) {
      j = i;
      while ((j < n*last) && (ReadEntry(file, value))) {
        val[j] = value;
        j+=last;
      }  

      // Exiting if the file is corrupted.

      if (j < (n*last)) {
        throw ArpackError(ArpackError::UNEXPECTED_EOF, "ARdfMatrix");
      }
    }

  }
  else {

    // Adjusting last if we are going to read the last block.

    if (last > n) {
      last  = n;
      block = 0;
    }
    last  -= first;
    strows = m;
    stcols = last;

    // Reading matrix data.

    j = 0;
    while ((j < m*last) && (ReadEntry(file, value))) {
      val[j++] = value;  
    }

    // Exiting if the file is corrupted.

    if (j < m*last) {
      throw ArpackError(ArpackError::UNEXPECTED_EOF, "ARdfMatrix");
    }
 
  }

} // ReadBlock.


template<class ARTYPE>
void ARdfMatrix<ARTYPE>::Define(const std::string& filename, int blksizep)
{

  // Declaring variables.

  char   c;
  char   data[LINELEN];

  // Opening the file.

  datafile = filename;
  file.open(datafile.c_str());
  
  if (!file) {
    throw ArpackError(ArpackError::CANNOT_OPEN_FILE, "ARdfMatrix");
  }

  // Setting initial values.

  blksize  = blksizep;
  block    = 0;
  headsize = 0;
  first    = 0;
  strows   = 0;
  stcols   = 0;

  // Reading the file heading.

  do {
    file.get((char*)data,LINELEN,'\n'); 
    file.get(c);
    headsize++;
  }
  while (data[0] == '%'); 

  // Reading m and n or returning if a problem was detected.

  if (sscanf(data, "%d %d", &m, &n) != 2) {
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARdfMatrix");
  }
  if ((m<1) || (n<1)) {
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARdfMatrix");
  }

  // Defining roword.

  roword = ((blksize != 0) && (m > n));

  // (Re)Dimensioning val.

  if (val != NULL) delete[] val;

  if (blksize == 0) {

    // Redefining blksize and reading the entire matrix.

    blksize = n;
    nblocks = 1;
    val = new ARTYPE[m*blksize]; 
    ReadBlock();
  
  }
  else if (roword) {  

    // m >> n, so we will read only blksize rows (but not now).
 
    if (blksize > m) blksize = m;
    nblocks = (m+blksize-1)/blksize;
    val = new ARTYPE[blksize*n]; 
    if (blksize == m) ReadBlock();

  } 
  else {       

    // n >> m, so we will read only blksize columns (but not now).
  
    if (blksize > n) blksize = n;
    nblocks = (n+blksize-1)/blksize;
    val = new ARTYPE[m*blksize]; 
    if (blksize == n) ReadBlock();

  }

} // Define.


template<class ARTYPE>
ARdfMatrix<ARTYPE>::ARdfMatrix()
{

  m        = 0; 
  n        = 0;
  block    = 0;
  blksize  = 0;
  headsize = 0;
  first    = 0;
  strows   = 0;
  stcols   = 0;
  roword   = false;
  val      = NULL;

} // Short constructor.


template<class ARTYPE>
ARdfMatrix<ARTYPE>::ARdfMatrix(const std::string& filename, int blksizep) 
{ 

  val = NULL;
  Define(filename, blksizep); 

} // Long constructor.


template<class ARTYPE>
ARdfMatrix<ARTYPE>::~ARdfMatrix()
{

  if (val != NULL) delete[] val;

} // Destructor.


#endif // ARDFMAT_H

