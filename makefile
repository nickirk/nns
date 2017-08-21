CXX=g++
LDFLAGS=-std=c++11 -g -Wall
EIGEN_PATH=lib/eigen/
EIGEN_FLAGS=-I$(EIGEN_PATH)
SOURCEFILES=Sampler.cxx Nnw.cxx Basis.cxx Hamiltonian.cxx Determinant.cxx
TESTFILES=testNnw.cxx testSampler.cxx testBasis.cxx
SRC=src
TST=test
BUILD=build
TSTBUILD=build/test
DDIR=.depend
TSTDDIR=.depend/test
SOURCES=$(patsubst %,$(SRC)/%,$(SOURCEFILES)) $(patsubst %,$(TST)/%,$(TESTFILES))
OBJECTS=$(patsubst %.cxx,$(BUILD)/%.o,$(SOURCEFILES)) 
DEPENDENCIES=$(patsubst %.cxx,$(DDIR)/%.d,$(SOURCEFILES)) $(patsubst %.cxx,$(TSTDDIR)/%.d,$(TESTFILES))

BTEST=basisTest
NNWTEST=nnwTest
STEST=samplerTest
HTEST=hamTest
ETEST=eigenTest
BTESTOBJECT=$(TSTBUILD)/testBasis.o
NNWTESTOBJECT=$(TSTBUILD)/testNnw.o
STESTOBJECT=$(TSTBUILD)/testSampler.o
ETESTOBJECT=$(TSTBUILD)/testEigen.o
HTESTOBJECT=$(TSTBUILD)/testHam.o

$(BUILD)/%.o: $(SRC)/%.cxx
	$(CXX) $(EIGEN_FLAGS) $(LDFLAGS) -c $< -o $@

$(DDIR)/%.d: $(SRC)/%.cxx
	$(CXX) $(EIGEN_FLAGS) $(LDFLAGS) -MM -MT $< -o $@

$(TSTBUILD)/%.o: $(TST)/%.cxx
	$(CXX) $(EIGEN_FLAGS) $(LDFLAGS) -c $< -o $@

$(TSTDDIR)/%.d: $(TST)/%.cxx
	$(CXX) $(EIGEN_FLAGS) $(LDFLAGS) -MM -MT $< -o $@

-include $(DEPENDENCIES)

testBasis: $(OBJECTS) $(BTESTOBJECT)
	$(CXX) $(EIGEN_FLAGS) $(LDLFLAGS) -o $(BTEST) $(OBJECTS) $(BTESTOBJECT)

testNnw: $(OBJECTS) $(NNWTESTOBJECT)
	$(CXX) $(EIGEN_FLAGS) $(LDLFLAGS) -o $(NNWTEST) $(OBJECTS) $(NNWTESTOBJECT)

testSampler: $(OBJECTS) $(STESTOBJECT)
	$(CXX) $(EIGEN_FLAGS) $(LDLFLAGS) -o $(STEST) $(OBJECTS) $(STESTOBJECT)

testEigen: $(OBJECTS) $(ETESTOBJECT)
	$(CXX) $(EIGEN_FLAGS) $(LDLFLAGS) -o $(ETEST) $(OBJECTS) $(ETESTOBJECT)

testHam: $(OBJECTS) $(HTESTOBJECT)
	$(CXX) $(EIGEN_FLAGS) $(LDLFLAGS) -o $(HTEST) $(OBJECTS) $(HTESTOBJECT)
