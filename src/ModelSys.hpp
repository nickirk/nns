#ifndef ModelSys_DEFINED
#define ModelSys_DEFINED

// type of the determiants
typedef std::vector<bool> detType;

inline int intcast(detType const & in);

class Basis{
  public:
    Basis(int numEle_, int numOrb_);
    int getSize();
  private:
    int numEle;
    int numOrb;
    int size;
    void generateBasis();
    void indexBasis();
    class BasisIndexRef();
};


class Hamiltonian{
  public:
    Hamiltonian(int size_, double diagDelta_, double offDiagMagntd_, double
		offDiagNonZeroRatio_){}
    int getSize() const;
    int getSparseSize() const;
    void sparseAcess(int pos, int &row, int &col, double &value);
    double & operator () (detType i, detType j);
    double & operator () (int i, int j);
    void multiplyMatrixVector(double *v, double *w);
  private:
    int size;
    double diagDelta;
    double offDiagMagntd;
    double OffDiagNonZeroRatio;
    std::vector<int> rowVec, colVec;
    std::vector<double> valueVec;
    void initHamiltonian();
};

#endif
