/*
 * AbInitioHamiltonian.h
 *
 *  Created on: Feb 13, 2018
 *      Author: Lauretta Schwarz
 */

#ifndef SRC_HAMILTONIAN_ABINITIOHAMILTONIAN_HPP_
#define SRC_HAMILTONIAN_ABINITIOHAMILTONIAN_HPP_

#include <string>
#include "Hamiltonian.hpp"

namespace networkVMC{

class excit_store;

class AbInitioHamiltonian: public Hamiltonian {
public:
    AbInitioHamiltonian(int dimension):Hamiltonian(dimension){psingles=0.0; pdoubles=0.0; pparallel=0.0;};
	virtual ~AbInitioHamiltonian();
	int getFermiSign(detType const &alpha, int annihilatorIndex, int creatorIndex) const;
    // count the number of connected states
    int countNumberCoupledStates(detType const &source, int exflag, int &nsingleexcit, int &ndoubleexcit);
    // calculate the generation probability
    double calcGenProp(detType const &source, detType const &target);
    // return the next single excitation (deterministic excitation generator)
    detType getSingleExcitation(detType const &source, std::vector<int> const &source_orbs, std::vector<int> &holes, std::vector<int> &particles, int &exflag, bool &ballexcitfound) const;
    // return the next double excitation (deterministic excitation generator)
    detType getDoubleExcitation(detType const &source, std::vector<int> const &source_orbs, std::vector<int> &holes, std::vector<int> &particles, int &exflag, bool &ballexcitfound) const;
    // randomly generate a single excitation excitation (random excitation generator)
    detType genSingleExcitation(detType const &source, std::vector<int> const &source_orbs, std::vector<int> &holes, std::vector<int> &particles, double &pGet, double pdoubnew, std::vector<int> noccs, std::vector<int> nunoccs) const;
    // randomly generate a double excitation excitation (random excitation generator)
    detType genDoubleExcitation(detType const &source, std::vector<int> const &source_orbs, std::vector<int> &holes, std::vector<int> &particles, double &pGet, double pdoubnew, std::vector<int> noccs, std::vector<int> nunoccs) const;
    // random excitation generator: randomly generates one connected determinant
    detType getRandomCoupledState(detType const &source, double &p) const;
    // deterministic excitation generator: generates all connected determinants
    std::vector<detType> getCoupledStates(detType const &source) const;
    // random excitation generator: randomly generates one connected determinant biased according to the 
    // connecting Hamiltonian matrix element
    detType getRandomCoupledState_cs(detType const &source, double &p) const;
    // calculate the generation probability biased according to the connecting
    // Hamiltonian matrix element
    double calcGenProp_cs(detType const &source, detType const &target) const;

    // a class for the excitation generation
    // forward declaration
    class excit_store;
    // probabilities for the random excitation generation
    // single and / or double excitations ?
    //int exflag;
    // probability of generating a single excitation
    double psingles;
    // probability of generating a double excitation
    double pdoubles;
    // probability of generating a parallel spin double excitation
    double pparallel;
    // dynamic adjustment of the probabilities
    bool bbias_sd;
    bool bbias_po;
    // which part of the excitation generator
    bool lin_exact = false;
    bool part_exact = true;
    // set values for the generation probabilities based on the number 
    // of electrons and holes
    void set_probabilities(detType example_det);
    // initialise these parameters
    void set_probabilities_bias(detType example_det, int exflag);
    // for updating these parameters
    void check_probabilities(excit_store const &excitation, double hel, int nspawns);
    //void check_probabilities(double hel, int nspawns);
    void adjust_probabilities();



private:
    // for the dynamic adjustment of the probabilities
    int minsingle;
    int mindouble;
    int minparadouble;
    int minoppdouble;
    double max_hel_single_ratio;
    double max_hel_double_ratio;
    double max_hel_para_double_ratio;
    double max_hel_opp_double_ratio;
    int nsingle;
    int ndouble;
    int nparadouble;
    int noppdouble;

    // for the excitation generation
    // return the next electron pair (for the deterministic excitation generator)
    void getElecPair(std::vector<int> const &source_orbs, int &el1, int &el2, int &spinpair, int &ind) const;
    // pick a random electron pair (random excitation generator)
    std::vector<int> pickElecPair(std::vector<int> const &source_orbs, std::vector<int> &spin, int &elecpairs) const;
    // could the number of forbidden spin orbitals (random excitation generator)
    int countForbiddenOrbs(std::vector<int> const &spins, std::vector<int> nunoccs) const;
    // pick orbital a of an ij->ab excitation (random excitation generator)
    int pickOrbA(detType const &source, std::vector<int> const &spins, std::vector<int> noccs, std::vector<int> unoccs, int norbs, int nel, int nforbiddenorbs, int &nexcita, int &aspin, bool &baorbfail) const;
    // pick orbital b of an ij->ab excitation (random excitation generator)
    int pickOrbB(detType const &source, std::vector<int> const &spins, std::vector<int> noccs, std::vector<int> unoccs, int norbs, int nel, int nforbiddenorbs, int aorb, int aspin, int &nexcitb, int &bspin, int &nexcitotherway) const;


};

// definition of nested class

class AbInitioHamiltonian::excit_store{

    public:

        // the source determinant in bit-string representation
        detType source;
        // source determinant in natural integer represenation
        std::vector<int> source_orbs;
        // excitation matrix
        int excitmat[2][2];
        // target determinant in bit-string representation
        detType target;
        // generation probability
        double pgen;
        // default constructor
        excit_store();
        // constructor
        excit_store(detType source_det){source=source_det; bfilled=false;};
        // fill in information necessary for random excitation generation
        void construct_class_count();
        // generate a single excitation
        detType gen_single_excit_cs(AbInitioHamiltonian const &H);
        // generate a double excitation
        detType gen_double_excit_cs(AbInitioHamiltonian const &H);
        // generate a random excitation
        void gen_rand_excit_cs(AbInitioHamiltonian const &H);
        // calculate the generation probability
        void calc_pgen_cs(AbInitioHamiltonian const &H, std::vector<int> const &src, std::vector<int> const &tgt);

    private:

        // indicate whether all variables are filled
        bool bfilled;
        // number of spin orbitals
        int norbs;
        // number of electrons
        int nel;
        // number of alpha spin electrons
        int nalphaels;
        // number of beta spin electrons
        int nbetaels;
        // number of alpha spin holes
        int nalphaholes;
        // number of beta spin holes
        int nbetaholes;
        // number of pairs of alpha-alpha spin electrons
        int aa_elec_pairs;
        // number of pairs of beta-beta spin electrons
        int bb_elec_pairs;
        // number of pairs of alpha-beta spin electrons
        int ab_elec_pairs;

        // select a hole for a single excitation
        int select_single_hole(AbInitioHamiltonian const &H, int src, double &pgen);
        // calculate the generation probability for a single excitation
        double pgen_single_excit_cs(AbInitioHamiltonian const &H, int src, int tgt);
        // calculate the generation probability for selecting holes in a doule excitation
        double pgen_select_double_hole(AbInitioHamiltonian const &H, std::vector<int> const &src, int const &orb_pair, double &cum_sum, int const &tgt);
        // pick a pair of electrons
        std::vector<int> pick_biased_elecs(AbInitioHamiltonian const &H, std::vector<int> &elecs);
        // pick a pair of holes (same spin pair)
        std::vector<int> select_double_holes(AbInitioHamiltonian const &H, std::vector<int> const &src, std::vector<double> &cum_sum, std::vector<double> &cpt);
        // pick a holes for an opposite spin pair
        int select_double_hole(AbInitioHamiltonian const &H, std::vector<int> const &src, int const &orb_pair, double &cum_sum, double &cpt);
        // contribution for a pair of opposite spin electrons
        double opp_spin_pair_contribution(AbInitioHamiltonian const &H, int i, int j, int a, int b);
        // contribution for a pair of same spin electrons
        double same_spin_pair_contribution(AbInitioHamiltonian const &H, int i, int j, int a, int b);

};

AbInitioHamiltonian readAbInitioHamiltonian(int dim, std::string file_name);

}

#endif /* SRC_HAMILTONIAN_ABINITIOHAMILTONIAN_HPP_ */
