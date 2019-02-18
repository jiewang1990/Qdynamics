//
//  quantumscar.h
//  quantumscar
//
//  Created by Jie Wang on 6/23/15.
//  Copyright (c) 2015 Jie Wang. All rights reserved.
//

#ifndef quantumscar_h
#define quantumscar_h

#include "phpfaffian.h"
////Spectra.
//#include <GenEigsSolver.h>
//#include <MatOp/SparseGenMatProd.h>
//#include <MatOp/SparseSymMatProd.h>
//#include <SymEigsSolver.h>

void quantum_scar_new(string filename);
void quantum_scar(string filename);
void quantum_scar_invt(string filename);
void quantum_scar_mps(string filename);
//void quantum_scar(int no, int np, int momentum, bool needM, bool useqh, int rangeS);
void plot_H_dim(int, int, string);

//Eigen::SparseMatrix<double> SparseRemoveRow(Eigen::SparseMatrix<double> input, int ind, vector<state_int>& bit_list, vector<state>& state_list);
Eigen::SparseMatrix<double> SparseLinearComb(Eigen::SparseMatrix<double> input, int ind1, int ind2, double coeff1, double coeff2);

struct eigen_set{
    bool parity;
    double particlehole;
    Eigen::VectorXd evec;
    double eval;
    double entanglement;
    double trans;
    complex<double> singletrans1, singletrans2;
};
inline bool compare_eigen_set(const eigen_set& inputa, const eigen_set& inputb) {
    double dif=inputa.eval-inputb.eval;
    if (abs(dif)>1e-12) {
        if (inputa.eval>inputb.eval) {
            return true;
        }
        else if (inputa.eval<inputb.eval) {
            return false;
        }
    }
    else {
        if (inputa.parity and !inputb.parity) {
            return true;
        }
        else if (!inputa.parity and inputb.parity) {
            return false;
        }
        else {
            return true;
        }
    }
}
inline bool compare_vecdouble(const vector<double>& a, const vector<double>& b) {
    if (a[0]<b[0]) {
        return true;
    }
    else return false;
}
inline char printparity(const bool& input) {
    if (input) return '+'; else return '-';
}
inline Eigen::MatrixXd commutator(Eigen::SparseMatrix<double> inputa, Eigen::SparseMatrix<double> inputb, bool commut){
    if (commut) {
        return inputa*inputb-inputb*inputa;
    }
    else {
        return Eigen::MatrixXd(inputa*inputb+inputb*inputa);
    }
}
//inline double safe_mult(double x, double y){
//    return x*y;
//}
//inline complex<double> safe_mult(complex<double> x, complex<double> y){
//    return x*conj(y);
//}

class Scars:public Msector{
private:
    int rangeS;
    int nvec1, nvec2;
    int replu_range, connectX;
    string bound_cond, type;
    bool fullsymmetry=true;//parity + T^(L/2).
    bool constrain;
    void init(int);
    void diag_setup(string);
    void make_mlist(string);
    void generate_hilberspace();
    vector<matrix> matrixlist_M_h1h1;
    inline bool legal_state(state_int);//to see if a state satisfies quantum-scar constrain: adjacent two sites has at least 0 for any site.
//    inline bool glegal_state(state_int);//generalized-legal_state, satisfies replusion constraint.
    //inline bool flip_state(state_int, int);//to see if to flip a spin or not.
    Eigen::SparseMatrix<double> H_spamatrix;
    
public:
    Scars();
    Scars(int, int, string, bool, string, int replu_range=1, int connectX=1);
    void out_info();
    void usefullsymmetry(bool);
    Msector Diag;
    void show_scar_energy(int);
    Eigen::SparseMatrix<double> get_H();
    vector<double> get_energies(string);
    void setup_symmetry();
    
    int get_H_dim();
    
    //assist.
    Eigen::SparseMatrix<double> SparseRemoveRow(Eigen::SparseMatrix<double> input, int ind, vector<state_int>& bit_list, vector<state>& state_list);
//    Eigen::SparseMatrix<double> SparseLinearComb(Eigen::SparseMatrix<double> input, int ind1, int ind2, double coeff1, double coeff2);
    
    //Inversion.
    void Inversion_setup(), proj_H_invsector(), show_scar_energy_inv(int, bool);
    void diag_scar_H();
    void diag_scar_H_inv(string, bool calculatetrans=false, bool calculatees=false, bool calculateph=false);//diag with parity.
    bool diag_scar_H_inv_halft(bool, bool, bool);//use parity and half-translation.
    bool diag_scar_H_inv_t(bool, int, bool calculateee=false);
    //void diag_scar_H_spectra(string, int, int);//diag with parity, and spectra.
    bool diag_scar(bool inv, int Ky, string mode, bool calculateee=false);
    vector<eigen_set> Eigen_Sets, Eigen_Sets_pls, Eigen_Sets_min;
    Eigen::SparseMatrix<double> Invmat;
    void make_Invmat();
    int Eval_N_pls, Eval_N_min;
    Eigen::MatrixXd P_pls, P_min, H_inv_pls, H_inv_min;
    
    //Eigen states.
    void Show_Eigen_Sets(int, bool printtrans=false, bool printes=false, bool printph=false);
    void Print_Eigen_Sets(int, string, bool printes=false, bool printph=false, bool printtrans=false);
    void Print_Eigen_Sets2(int, string, bool printes=false, bool printtrans=false);
    void pick_Eigen_Sets(double E1, double E2, double S1, double S2, vector<eigen_set>& eigenset, vector<Eigen::VectorXd>& spectrum);
    
    //project into Inv+Trans sector.
    Eigen::MatrixXd reducedH;
    bool make_invtrans_proj(bool, bool, bool zerosector=false);
    //bool make_invtrans_proj_0sector(bool);
    Eigen::SparseMatrix<double> invtrans_proj;
    bool make_inv_trans0pi_proj(bool, bool);
    Eigen::MatrixXd inv_trans0pi_proj;
    //vector<state> statelist_0pi;
    void make_inversion_proj();
    Eigen::SparseMatrix<double> shrinkMatrix_trans;
    Eigen::SparseMatrix<double> shrinkMatrix_trans_inv;
    
    //make scar shrinker.
    void makeScarShrinker_trans(int);//project into 0 or pi sector.
    void makeScarShrinker_trans_inv(bool);
    vector<state> statelist_M_trans_int;
    vector<state_int> bitlist_trans_int;
    void print_bitlist_trans_int();
    
    //Particle-Hole.
    void ParticleHole_setup();
    Eigen::SparseMatrix<double> PHmat;
    
    //Translation.
    void Translation_setup();
    Eigen::SparseMatrix<double> Transmat;
    Eigen::SparseMatrix<double> Transmat_halfL;//=Trans_pow(L/2) if L even; =Trans_pow(L) if L odd.
    Eigen::SparseMatrix<double> Trans_pow(int);
    
    void print_commutators();
    
    //Time evolution.
    void Time_Evolution(const double dt, const int Nt, const state_int& state);
    
    //MPS test state.
    Eigen::MatrixXd B0 ,B1, C0, C1;
    void mps_run(int);
    void mps_run();
    void test_mps_wfs();
    Eigen::VectorXd mpswf;
    Eigen::VectorXd mpswf0, mpswf1, mpswf_pls, mpswf_min;
    void print_mpswf();
    vector<Eigen::VectorXd> mps_entanglement();
    
    ~Scars();
    
};

#endif
