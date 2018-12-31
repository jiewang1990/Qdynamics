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

void quantum_scar(string filename);
//void quantum_scar(int no, int np, int momentum, bool needM, bool useqh, int rangeS);

struct eigen_set{
    bool parity;
    double particlehole;
    Eigen::VectorXd evec;
    double eval;
    double entanglement;
    double trans;
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
    string bound_cond, type;
    bool fullsymmetry=true;
    void init(int);
    void diag_setup(string);
    void make_mlist(string);
    void make_mlist_PXP();
    void make_mlist_PX2P();
    void generate_scar_hilberspace();
    vector<matrix> matrixlist_M_h1h1;
    inline bool legal_state(state_int);//to see if a state satisfies quantum-scar constrain: adjacent two sites has at least 0 for any site.
    //inline bool flip_state(state_int, int);//to see if to flip a spin or not.
    Eigen::SparseMatrix<double> H_spamatrix;
    
public:
    Scars();
    Scars(int, int, string, string, int, int);
    void usefullsymmetry(bool);
    Msector Diag;
    void diag_scar();
    void show_scar_energy(int);
    Eigen::SparseMatrix<double> get_H();
    vector<double> get_energies(string);
    void setup_symmetry();
    
    //Inversion.
    void Inversion_setup(), proj_H_invsector(), show_scar_energy_inv(int, bool);
    void Show_Eigen_Sets(int, bool printtrans=false, bool printes=false, bool printph=false);
//    void Show_Eigen_Sets_ph(int, bool printes=false);
    void Print_Eigen_Sets(int, string, bool printes=false, bool printph=false);
    void diag_scar_H_inv(string, bool calculatetrans=false, bool calculatees=false, bool calculateph=false);//diag with parity.
    bool diag_fullsymmetry(bool, bool, bool);//use parity and translation.
//    void diag_scar_H_spectra(string, int, int);//diag with parity, and spectra.
    vector<eigen_set> Eigen_Sets, Eigen_Sets_pls, Eigen_Sets_min;
    Eigen::SparseMatrix<double> Invmat;
    void make_Invmat();
    int Eval_N_pls, Eval_N_min;
    Eigen::MatrixXd P_pls, P_min, H_inv_pls, H_inv_min;
    
    //project into Inv+Trans sector.
    void proj_full_sym();
    Eigen::MatrixXd reducedH;
    bool make_invtrans_proj(bool, bool);
    Eigen::MatrixXd invtrans_proj;
    void make_inversion_proj();
    
    //Particle-Hole.
    void ParticleHole_setup();
    Eigen::SparseMatrix<double> PHmat;
    
    //Translation.
    void Translation_setup();
    Eigen::SparseMatrix<double> Transmat;
    Eigen::SparseMatrix<double> Transmat_halfL;//=Trans_pow(L/2) if L even; =Trans_pow(L) if L odd.
    Eigen::SparseMatrix<double> Trans_pow(int);
    
    void print_commutators();
    
    //MPS test state.
    Eigen::MatrixXd B0 ,B1, C0, C1;
    void mps_run(int);
    void test_mps_wfs();
    Eigen::VectorXd mpswf;
    void print_mpswf();
    
    ~Scars();
    
};

#endif
