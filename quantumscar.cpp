//
//  quantumscar.cpp
//  quantumscar
//
//  Created by Jie Wang on 12/10/18.
//  Copyright © 2018 Jie Wang. All rights reserved.
//

#include "quantumscar.h"

Scars::Scars(){};
Scars::Scars(int no_t, int ranges_t, string bound_cond_t, int nvec1_t, int nvec2_t):rangeS(ranges_t), bound_cond(bound_cond_t), nvec1(nvec1_t), nvec2(nvec2_t) {
    this->init(no_t);
    this->generate_scar_hilberspace();
    this->setup_symmetry();
    
    this->diag_setup();
    
    cout<<"$ project H into inverse sector"<<endl;
    this->proj_H_invsector();

//    cout<<"$ start diagonalization $"<<endl;
//    //this->diag_scar();
//    this->diag_scar_H_inv("both");
////    this->diag_scar_H_spectra("both", this->nvec1, this->nvec2);
//    cout<<"$ done  diagonalization $"<<endl;
//
//    this->print_Eigen_Sets(200);
}
void Scars::init(int no_t){
    this->statelist_M.clear(); this->matrixlist_M.clear(); this->one=(state_int) 1;
    //TODO: only for spin-half, it is "fermion".
    this->No=no_t; this->Np=0; this->momentum=0; this->geometry="sphere"; this->statistics="fermion"; this->shift=vector<double>(2, 0.); this->useqh=false; this->needM=false;
}
void Scars::setup_symmetry(){
    if (this->bound_cond=="PBC") {
        this->Inversion_setup();
        this->ParticleHole_setup();
        this->Translation_setup();
    }
    else if (this->bound_cond=="OBC") {
        this->Inversion_setup();
        this->ParticleHole_setup();
    }
}
void Scars::diag_setup(){
    this->make_mlist_h1h1();
    this->matrixlist_M=mergemlist(sortmlist(this->matrixlist_M));
    this->H_spamatrix=Eigen::SparseMatrix<double>(this->bitlist.size(), this->bitlist.size());
    vector<Eigen::Triplet<double>> triplets=mlist_to_ftriplet(this->matrixlist_M);
    this->H_spamatrix.setFromTriplets(triplets.begin(), triplets.end());
}
void Scars::diag_setup_perturb(){
    //TODO: make h1 diag. stored in mlist_h1.
    
}
inline bool Scars::legal_state(state_int input){
    if (bound_cond=="PBC") {
        for (int i=0; i<this->No; i++) {
            if (this->one<<supermod(i, this->No)&input and this->one<<supermod(i+1, this->No)&input) {
                return false;
            }
        }
        return true;
    }
    else if (bound_cond=="OBC") {
        //(i,i+i) w/ i<=No-2.
        for (int i=0; i<this->No-1; i++) {
            if (this->one<<supermod(i, this->No)&input and this->one<<supermod(i+1, this->No)&input) {
                return false;
            }
        }
        return true;
    }
    else {
        cout<<"unrecognized boundary condition"<<endl;
        exit(0);
    }
}
//inline bool Scars::flip_state(state_int input, int i){
//    if (this->bound_cond=="PBC") {
//        if (this->one<<supermod(i-1, this->No)&input and this->one<<supermod(i+1, this->No)&input) {
//            return false;
//        }
//        return true;
//    }
//    else if (this->bound_cond=="OBC") {
//        if (i!=0 and i!=this->No-1) {
//            if (this->one<<supermod(i-1, this->No)&input and this->one<<supermod(i+1, this->No)&input) {
//                return false;
//            }
//            return true;
//        }
//        else if (i==0) {
//            if (this->one<<supermod(i+1, this->No)&input) {
//                return false;
//            }
//            return true;
//        }
//        else if (i==this->No-1) {
//            if (this->one<<supermod(i-1, this->No)&input) {
//                return false;
//            }
//            return true;
//        }
//    }
//    else {
//        cout<<"unrecognized boundary condition"<<endl;
//        exit(0);
//    }
//}
void Scars::generate_scar_hilberspace(){
    this->basis.clear();
    for (state_int i=0; i<this->one<<this->No; i++) {
        if (this->legal_state(i)) {
            this->basis.push_back(i);
        }
    }
    //this->sort_f();
    this->statelist_M=vector<state>(this->basis.size());
    for (unsigned int i=0; i<this->basis.size(); i++) {
        this->statelist_M[i].obasis=int_to_vec(this->basis[i], this->No);
        this->statelist_M[i].Mom=0;
    }
    this->bitlist=this->basis;
}
void Scars::make_mlist_h0(){
    this->matrixlist_M.clear();
    
    for (unsigned int i=0; i<this->bitlist.size(); i++) {
        matrix k;
        k.ket=i;
        state_int bit_tmp0=this->bitlist[i], bit_tmp;
        for (int j=0; j<this->No; j++) {
            
            bit_tmp=bit_tmp0 ^ this->one<<j;
            if (!this->legal_state(bit_tmp)) {
                //print_bit(bit_tmp, this->No); cout<<"not legal state"<<endl;
                continue;
            }
            
            vector<state_int>::iterator binaryit=this->bitlist.begin();
            binaryit=lower_bound(this->bitlist.begin(), this->bitlist.end(), bit_tmp);
            if (binaryit!=this->bitlist.end() and *binaryit==bit_tmp) {
                this->index=binaryit-this->bitlist.begin();
                k.bra=index;
                k.element=1.0;
                this->matrixlist_M.push_back(k);
                //print_bit(this->bitlist[i], this->No); print_bit(bit_tmp, this->No); cout<<"ind="<<index<<endl<<endl;
            }
            else {
                cout<<"cannot find basis"<<endl;
                print_bit(this->bitlist[i], this->No); print_bit(bit_tmp, this->No);
                exit(0);
            }
        }
    }
}
void Scars::make_mlist_h1h1(){
    this->matrixlist_M.clear();//TODO: to change !!!!!!
    
    for (unsigned i=0; i<this->bitlist.size(); i++) {
        matrix k;
        k.ket=i;
        state_int bit_tmp0=this->bitlist[i], bit_tmp;
        for (int j1=0; j1<this->No; j1++) {
            for (int j2=j1; j2<this->No; j2++) {
                
                bit_tmp=bit_tmp0 ^ this->one<<j1;
                bit_tmp=bit_tmp  ^ this->one<<j2;
                
                if (!this->legal_state(bit_tmp)) {
                    continue;
                }
                
                vector<state_int>::iterator binaryit=this->bitlist.begin();
                binaryit=lower_bound(this->bitlist.begin(), this->bitlist.end(), bit_tmp);
                if (binaryit!=this->bitlist.end() and *binaryit==bit_tmp) {
                    this->index=binaryit-this->bitlist.begin();
                    k.bra=index;
                    k.element=2.0;
                    if (j1==j2) {
                        k.element=1.0;
                    }
                    this->matrixlist_M.push_back(k);
                    //print_bit(this->bitlist[i], this->No); print_bit(bit_tmp, this->No); cout<<"ind="<<index<<endl<<endl;
                }
                else {
                    cout<<"cannot find basis"<<endl;
                    print_bit(this->bitlist[i], this->No); print_bit(bit_tmp, this->No);
                    exit(0);
                }
            }
        }
    }
}
void Scars::diag_scar(){
    cout<<"@@ begin diagonalize @@"<<endl;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,-1,-1>> es;
    es.compute(this->H_spamatrix);
    this->Eval=vector<double>(es.eigenvalues().size(), 0.);
    for (int i=0; i<Eval.size(); i++) {
        this->Eval[i]=es.eigenvalues()[i];
    }
    this->Ground_State=Eigen::VectorXcd::Zero(this->statelist_M.size());
    for (int i=0; i<this->Ground_State.size(); i++) {
        this->Ground_State(i)=es.eigenvectors().col(0)[i];
    }
    cout<<"@@ done  diagonalize @@"<<endl;
}
Eigen::SparseMatrix<double> Scars::get_H(){
    return this->H_spamatrix;
}
void Scars::show_scar_energy(int range){
    int N=this->Eval.size();
    range/=2;
    cout<<"\n&&& energies &&&"<<endl;
    if (N%2==1) {
        N=(N-1)/2;
        range=min(range, N);
        for (int i=N-range; i<=N+range; i++) {
            cout<<setprecision(6)<<chop(this->Eval[i])<<endl;
        }
    }
    else {
        N=N/2;
        range=min(range-1, N-1);
        for (int i=N-range-1; i<=N+range; i++) {
            cout<<setprecision(6)<<chop(this->Eval[i])<<endl;
        }
    }
    cout<<"&&&   done   &&&\n"<<endl;
    //The dimension of zero modes is: L%2==1, Fib[(L-1)/2]; L%2==0, Fib[L/2+1].
}

//Inversion.
void Scars::Inversion_setup(){
    this->Invmat=Eigen::SparseMatrix<double>(this->bitlist.size(), this->bitlist.size());
    
    vector<Eigen::Triplet<double>> triples; triples.clear();
    
    for (unsigned int i=0; i<this->bitlist.size(); i++) {
        state_int bit_tmp=invert_bits(this->bitlist[i], this->No);
        vector<state_int>::iterator binaryit=this->bitlist.begin();
        binaryit=lower_bound(this->bitlist.begin(), this->bitlist.end(), bit_tmp);
        
        if (binaryit!=this->bitlist.end() and *binaryit==bit_tmp) {
            triples.push_back(Eigen::Triplet<double>(i, binaryit-this->bitlist.begin(), 1));
            //print_bit(this->bitlist[i], this->No); print_bit(bit_tmp, this->No); cout<<"ind="<<index<<endl<<endl;
        }
        else {
            cout<<"cannot find inverted bit"<<endl;
            exit(0);
        }
        
    }
    
    this->Invmat.setFromTriplets(triples.begin(), triples.end());
    
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,-1,-1>> es;
    es.compute(this->Invmat);
    vector<double> InvEval=vector<double>(es.eigenvalues().size(), 0.);
    for (int i=0; i<InvEval.size(); i++) {
        InvEval[i]=es.eigenvalues()[i];
    }
    vector<Eigen::VectorXd> InvEvec=vector<Eigen::VectorXd>(es.eigenvalues().size());
    for (int i=0; i<InvEval.size(); i++) {
        InvEvec[i]=Eigen::VectorXd::Zero(this->statelist_M.size());
        for (int j=0; j<this->bitlist.size(); j++) {
            InvEvec[i](j)=es.eigenvectors().col(i)[j];
        }
    }
    
    this->Eval_N_pls=0; this->Eval_N_min=0;
    for (int i=0; i<InvEval.size(); i++) {
        if (InvEval[i]>0) {
            Eval_N_pls++;
        }
        else {
            Eval_N_min++;
        }
    }
    
    this->P_pls=Eigen::MatrixXd::Zero(this->bitlist.size(), Eval_N_pls);
    this->P_min=Eigen::MatrixXd::Zero(this->bitlist.size(), Eval_N_min);
    int counter_pls=0, counter_min=0;
    for (int i=0; i<InvEval.size(); i++) {
        if (InvEval[i]>0) {
            this->P_pls.col(counter_pls++)=InvEvec[i];
        }
        else if (InvEval[i]<0) {
            this->P_min.col(counter_min++)=InvEvec[i];
        }
    }
}
void Scars::proj_H_invsector(){
    this->H_inv_pls=this->P_pls.adjoint()*this->H_spamatrix*this->P_pls;
    this->H_inv_min=this->P_min.adjoint()*this->H_spamatrix*this->P_min;
}
void Scars::diag_scar_H_inv(string invpls, bool calculatees, bool calculateph){
    if (invpls=="pls") {
        if (this->H_inv_pls.size()==0) {
            cout<<"to make H_inv_pls"<<endl;
            exit(0);
        }
        else {
            //cout<<"@@ begin diagonalize (+) @@"<<endl;
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,-1,-1>> es;
            es.compute(this->H_inv_pls);
            
            this->Eigen_Sets_pls.clear();
            for (int i=0; i<es.eigenvalues().size(); i++) {
                eigen_set eset_tmp;
                eset_tmp.parity=true;
                eset_tmp.evec=es.eigenvectors().col(i);
                eset_tmp.eval=es.eigenvalues()[i];
                if (calculatees) {
                    Eigen::MatrixXd rho2;
                    this->ee_compute_rho(Xcd_to_Xd(this->P_pls*Xd_to_Xcd(eset_tmp.evec)), rho2, this->bitlist);
                    eset_tmp.entanglement=this->ee_eval_rho(rho2);
                }
                if (calculateph) {
                    Eigen::VectorXd vectmp=Xcd_to_Xd(this->P_pls*Xd_to_Xcd(eset_tmp.evec));
                    eset_tmp.particlehole=vectmp.dot(this->PHmat*vectmp);
                }
                this->Eigen_Sets_pls.push_back(eset_tmp);
            }
        }
    }
    else if (invpls=="min") {
        if (this->H_inv_min.size()==0) {
            cout<<"to make H_inv_min"<<endl;
            exit(0);
        }
        else {
            //cout<<"@@ begin diagonalize (-) @@"<<endl;
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,-1,-1>> es;
            es.compute(this->H_inv_min);
            
            this->Eigen_Sets_min.clear();
            for (int i=0; i<es.eigenvalues().size(); i++) {
                eigen_set eset_tmp;
                eset_tmp.parity=false;
                eset_tmp.evec=es.eigenvectors().col(i);
                eset_tmp.eval=es.eigenvalues()[i];
                if (calculatees) {
                    Eigen::MatrixXd rho2;
                    this->ee_compute_rho(Xcd_to_Xd(this->P_min*Xd_to_Xcd(eset_tmp.evec)), rho2, this->bitlist);
                    eset_tmp.entanglement=this->ee_eval_rho(rho2);
                }
                if (calculateph) {
                    Eigen::VectorXd vectmp=Xcd_to_Xd(this->P_min*Xd_to_Xcd(eset_tmp.evec));
                    eset_tmp.particlehole=vectmp.dot(this->PHmat*vectmp);
                }
                this->Eigen_Sets_min.push_back(eset_tmp);
            }
        }
    }
    else if (invpls=="both") {
        this->diag_scar_H_inv("pls", calculatees, calculateph);
        this->diag_scar_H_inv("min", calculatees, calculateph);
        this->Eigen_Sets.clear();
        this->Eigen_Sets.insert(this->Eigen_Sets.end(), this->Eigen_Sets_pls.begin(), this->Eigen_Sets_pls.end());
        this->Eigen_Sets.insert(this->Eigen_Sets.end(), this->Eigen_Sets_min.begin(), this->Eigen_Sets_min.end());
        
        sort(this->Eigen_Sets.begin(), this->Eigen_Sets.end(), compare_eigen_set);
    }
}
//void Scars::diag_scar_H_spectra(string invpls, int nvec1, int nvec2){
//    if (invpls=="pls") {
//        if (this->H_inv_pls.size()==0) {
//            cout<<"to make H_inv_pls"<<endl;
//            exit(0);
//        }
//        else {
//            Eigen::SparseMatrix<double> spamatrix=this->H_inv_pls.sparseView();
//            Spectra::SparseSymMatProd<double> op(spamatrix);
//            Spectra::SymEigsSolver<double, Spectra::SMALLEST_MAGN, Spectra::SparseSymMatProd<double>> es(&op, nvec1, nvec2);
//            cout<<"to inialization"<<endl;
//            es.init();
//            cout<<"to compute"<<endl;
//            int nconv = es.compute();
//            cout<<"spamatrix.cols()="<<spamatrix.cols()<<"   spamatrix.rows()="<<spamatrix.rows()<<endl;
//            cout<<"nonzeros="<<spamatrix.nonZeros()<<endl;
//            if (es.info() == Spectra::SUCCESSFUL) {
//                this->Eigen_Sets_pls.clear();
//                for (int i=0; i<es.eigenvalues().size(); i++) {
//                    eigen_set eset_tmp;
//                    eset_tmp.parity=true;
//                    eset_tmp.evec=es.eigenvectors().col(i);
//                    eset_tmp.eval=es.eigenvalues()[i];
//                    this->Eigen_Sets_pls.push_back(eset_tmp);
//                }
//            }
//            else {
//                cout<<"Not Successful."<<endl;
//                exit(0);
//            }
//        }
//    }
//    else if (invpls=="min") {
//        if (this->H_inv_min.size()==0) {
//            cout<<"to make H_inv_min"<<endl;
//            exit(0);
//        }
//        else {
//            Eigen::SparseMatrix<double> spamatrix=this->H_inv_min.sparseView();
//            Spectra::SparseSymMatProd<double> op(spamatrix);
//            Spectra::SymEigsSolver<double, Spectra::SMALLEST_MAGN, Spectra::SparseSymMatProd<double>> es(&op, nvec1, nvec2);
//            es.init();
//            int nconv = es.compute();
//            cout<<"spamatrix.cols()="<<spamatrix.cols()<<"   spamatrix.rows()="<<spamatrix.rows()<<endl;
//            if (es.info() == Spectra::SUCCESSFUL) {
//                this->Eigen_Sets_min.clear();
//                for (int i=0; i<es.eigenvalues().size(); i++) {
//                    eigen_set eset_tmp;
//                    eset_tmp.parity=false;
//                    eset_tmp.evec=es.eigenvectors().col(i);
//                    eset_tmp.eval=es.eigenvalues()[i];
//                    this->Eigen_Sets_min.push_back(eset_tmp);
//                }
//            }
//            else {
//                cout<<"Not Successful."<<endl;
//                exit(0);
//            }
//        }
//    }
//    else if (invpls=="both") {
//        this->diag_scar_H_spectra("pls", nvec1, nvec2);
//        this->diag_scar_H_spectra("min", nvec1, nvec2);
//        this->Eigen_Sets.clear();
//        this->Eigen_Sets.insert(this->Eigen_Sets.end(), this->Eigen_Sets_pls.begin(), this->Eigen_Sets_pls.end());
//        this->Eigen_Sets.insert(this->Eigen_Sets.end(), this->Eigen_Sets_min.begin(), this->Eigen_Sets_min.end());
//
//        sort(this->Eigen_Sets.begin(), this->Eigen_Sets.end(), compare_eigen_set);
//    }
//}
void Scars::Show_Eigen_Sets(int range, bool printes, bool printph){
    range/=2;
    int N=this->Eigen_Sets.size();
    if (N%2==1) {
        N=(N-1)/2;
        range=min(range, N);
        for (int i=N-range; i<=N+range; i++) {
            cout<<setprecision(6)<<setw(15)<<chop(this->Eigen_Sets[i].eval)<<",   "<<printparity(this->Eigen_Sets[i].parity);
            if (printph) {
                cout<<",   "<<setw(6)<<chop(this->Eigen_Sets[i].particlehole);
            }
            if (printes) {
                cout<<",   "<<setw(15)<<chop(this->Eigen_Sets[i].entanglement);
            }
            cout<<endl;
        }
    }
    else {
        N=N/2;
        range=min(range-1, N-1);
        for (int i=N-range-1; i<=N+range; i++) {
            cout<<setprecision(6)<<setw(15)<<chop(this->Eigen_Sets[i].eval)<<",   "<<printparity(this->Eigen_Sets[i].parity);
            if (printph) {
                cout<<",   "<<setw(6)<<chop(this->Eigen_Sets[i].particlehole);
            }
            if (printes) {
                cout<<",   "<<setw(15)<<chop(this->Eigen_Sets[i].entanglement);
            }
            cout<<endl;
        }
    }
    
    int zeromodedim=0;
    for (int i=0; i<this->Eigen_Sets.size(); i++) {
        if (abs(this->Eigen_Sets[i].eval)<1e-12) {
            zeromodedim++;
        }
    }
    cout<<"Zero mode dimension = "<<zeromodedim<<endl;
}
void Scars::Print_Eigen_Sets(int range, string filename, bool printes, bool printph){
    ofstream outfile(filename);
    
    range/=2;
    int N=this->Eigen_Sets.size();
    if (N%2==1) {
        N=(N-1)/2;
        range=min(range, N);
        for (int i=N-range; i<=N+range; i++) {
            outfile<<setprecision(6)<<setw(15)<<chop(this->Eigen_Sets[i].eval)<<"   "<<this->Eigen_Sets[i].parity;
            if (printph) {
                outfile<<"   "<<setw(6)<<chop(this->Eigen_Sets[i].particlehole);
            }
            if (printes) {
                outfile<<"   "<<setw(15)<<chop(this->Eigen_Sets[i].entanglement);
            }
            outfile<<endl;
        }
    }
    else {
        N=N/2;
        range=min(range-1, N-1);
        for (int i=N-range-1; i<=N+range; i++) {
            outfile<<setprecision(6)<<setw(15)<<chop(this->Eigen_Sets[i].eval)<<"   "<<this->Eigen_Sets[i].parity;
            if (printph) {
                outfile<<"   "<<setw(6)<<chop(this->Eigen_Sets[i].particlehole);
            }
            if (printes) {
                outfile<<"   "<<setw(15)<<chop(this->Eigen_Sets[i].entanglement);
            }
            outfile<<endl;
        }
    }
    outfile.close();
}

//Particle-Hole.
void Scars::ParticleHole_setup(){
    this->PHmat=Eigen::SparseMatrix<double>(this->bitlist.size(), this->bitlist.size());
    vector<Eigen::Triplet<double>> triples; triples.clear();
    
    for (unsigned int i=0; i<this->bitlist.size(); i++) {
        double sign=1.0;
        if ((this->No-count_bits(this->bitlist[i]))%2==1) {
            sign=-1.0;
        }
        triples.push_back(Eigen::Triplet<double>(i, i, sign));
    }
    
    this->PHmat.setFromTriplets(triples.begin(), triples.end());
}
//Translation.
void Scars::Translation_setup(){
    this->Transmat=Eigen::SparseMatrix<double>(this->bitlist.size(), this->bitlist.size());
    vector<Eigen::Triplet<double>> triples; triples.clear();
    vector<state_int>::iterator binaryit=this->bitlist.begin();
    
    for (unsigned int i=0; i<this->bitlist.size(); i++) {
        state_int statein=this->bitlist[i], stateout=0; int sign=1;
        
        cycle_bit(statein, stateout, this->No, 1, sign);
        
        binaryit=lower_bound(this->bitlist.begin(), this->bitlist.end(), stateout);
        
        if (binaryit!=this->bitlist.end() and *binaryit==stateout) {
            triples.push_back(Eigen::Triplet<double>(i, binaryit-this->bitlist.begin(), 1));
            //print_bit(this->bitlist[i], this->No); print_bit(stateout, this->No); cout<<"ind="<<binaryit-this->bitlist.begin()<<endl<<endl;
        }
        else {
            print_bit(this->bitlist[i], this->No); print_bit(stateout, this->No); cout<<"ind="<<binaryit-this->bitlist.begin()<<endl<<endl;
            cout<<"cannot find translated bit"<<endl;
            exit(0);
        }
    }
    
    this->Transmat.setFromTriplets(triples.begin(), triples.end());
    //cout<<"Transmat=\n"<<this->Transmat<<endl;
}
Eigen::SparseMatrix<double> Scars::Trans_pow(int x){
    if (x<1) {
        cout<<"cannot run Trans_pow"<<endl;
        exit(0);
    }
    else {
        Eigen::SparseMatrix<double> output=this->Transmat;
        for (int i=1; i<x; i++) {
            output=output*this->Transmat;
        }
        return output;
    }
}
void Scars::print_commutators(){
    cout<<"[H, P]=\n"<<commutator(this->H_spamatrix, this->Invmat, true)<<endl;
    cout<<"[H, C]=\n"<<commutator(this->H_spamatrix, this->PHmat,  true)<<endl;
    cout<<"{H, C}=\n"<<commutator(this->H_spamatrix, this->PHmat, false)<<endl;
    cout<<"[C, P]=\n"<<commutator(this->PHmat, this->Invmat, true)<<endl;
    cout<<"{C, P}=\n"<<commutator(this->PHmat, this->Invmat, false)<<endl;
    
    if (this->bound_cond=="PBC") {
        cout<<"[P, T]=\n"<<commutator(this->Invmat, this->Transmat, true)<<endl;
        cout<<"{P, T}=\n"<<commutator(this->Invmat, this->Transmat, false)<<endl;
        
        cout<<"[P, T^L]=\n"<<commutator(this->Invmat, this->Trans_pow(this->No), true)<<endl;
        //cout<<"{P, T^L}=\n"<<commutator(this->Invmat, this->Trans_pow(this->No), false)<<endl;
        
        cout<<"[P, T^(L/2)]=\n"<<commutator(this->Invmat, this->Trans_pow(this->No/2), true)<<endl;
        //cout<<"{P, T^(L/2)}=\n"<<commutator(this->Invmat, this->Trans_pow(this->No/2), false)<<endl;
        
        cout<<"[C, T]=\n"<<commutator(this->PHmat, this->Transmat, true)<<endl;
        cout<<"{C, T}=\n"<<commutator(this->PHmat, this->Transmat, false)<<endl;
        
        cout<<"[C, T^L]=\n"<<commutator(this->PHmat, this->Trans_pow(this->No), true)<<endl;
        cout<<"{C, T^L}=\n"<<commutator(this->PHmat, this->Trans_pow(this->No), false)<<endl;
        
        cout<<"[C, T^(L/2)]=\n"<<commutator(this->PHmat, this->Trans_pow(this->No/2), true)<<endl;
        cout<<"{C, T^(L/2)}=\n"<<commutator(this->PHmat, this->Trans_pow(this->No/2), false)<<endl;
    }
}
//MPS test state.
void Scars::mps_run(int ind){
    if (this->No%2==1) {
        cout<<"L%2==1, quit"<<endl;
        exit(0);
    }
    
    if (this->bound_cond=="OBC") {
        B0=Eigen::MatrixXd(2, 3);
        B1=Eigen::MatrixXd(2, 3);
        C0=Eigen::MatrixXd(3, 2);
        C1=Eigen::MatrixXd(3, 2);
        Eigen::VectorXd v=Eigen::VectorXd(2);
        if (ind==0) {
            v<<1., 1.;
        }
        else if (ind==1) {
            v<<1., -1.;
        }
        else if (ind==2) {
            v<<1., 0.;
        }
        else if (ind==3) {
            v<<0., 1.;
        }
        else {
            cout<<"cannot setup v"<<endl;
            exit(0);
        }
        
        B0  <<  1,  0,  0,  0,       1,  0;
        B1  <<  0,  0,  0,  sqrt(2.),0,  sqrt(2.);
        C0  <<  0, -1,  1,  0,       0,  0;
        C1  <<sqrt(2.),0,0, 0,-sqrt(2.), 0;
        
        //cout<<"B0=\n"<<B0<<endl;cout<<"B1=\n"<<B1<<endl;
        //cout<<"C0=\n"<<C0<<endl;cout<<"C1=\n"<<C1<<endl;
        //cout<<"v=\n"<<v<<endl;
        
        this->mpswf=Eigen::VectorXd::Zero(this->bitlist.size());
//        Eigen::MatrixXd matrixprod=v.transpose();
//        Eigen::VectorXd matrixprod=v.transpose();
        
//        cout<<"matrixprod=\n"<<matrixprod<<endl;
//        matrixprod*=B0;
//        cout<<"matrixprod=\n"<<matrixprod<<endl;
//        matrixprod*=C0;
//        cout<<"matrixprod=\n"<<matrixprod<<endl;
        
//        cout<<endl<<endl;
        
        for (int i=0; i<this->bitlist.size(); i++) {
            Eigen::MatrixXd matrixprod=v.transpose();
            for (int j=0; j<this->No/2; j++) {
                if (this->statelist_M[i].obasis[2*j]==0) {
                    matrixprod*=B0;
                }
                else {
                    matrixprod*=B1;
                }
                if (this->statelist_M[i].obasis[2*j+1]==0) {
                    matrixprod*=C0;
                }
                else {
                    matrixprod*=C1;
                }
            }
            matrixprod*=v;
            if (matrixprod.size()!=1) {
                cout<<"matrixprod.size()!=1"<<endl;
                exit(0);
            }
            this->mpswf(i)=matrixprod(0, 0);
        }
    }
    else {
        cout<<"cannot do PBC so far"<<endl;
        exit(0);
    }
    
    this->mpswf/=this->mpswf.norm();
    
    cout<<"$$ test the <H> and <H2> $$"<<endl;
    cout<<"<H> ="<<chop(this->mpswf.transpose()*this->H_spamatrix*this->mpswf)<<endl;
    cout<<"<H2>="<<chop(this->mpswf.transpose()*this->H_spamatrix*this->H_spamatrix*this->mpswf)<<endl;
    cout<<"$$ done the <H> and <H2> $$"<<endl;
}
void Scars::test_mps_wfs(){
    this->mps_run(0);
    Eigen::VectorXd mpswf1=this->mpswf;
    cout<<"<P> ="<<this->mpswf.dot(this->Invmat*this->mpswf)<<endl;
    cout<<"<P2>="<<this->mpswf.dot(this->Invmat*this->Invmat*this->mpswf)<<endl;
    cout<<"<PH> ="<<this->mpswf.dot(this->PHmat*this->mpswf)<<endl;
    cout<<"<PH2>="<<this->mpswf.dot(this->PHmat*this->Invmat*this->mpswf)<<endl;
    this->print_mpswf();
    
    this->mps_run(1);
    Eigen::VectorXd mpswf2=this->mpswf;
    cout<<"<P> ="<<this->mpswf.dot(this->Invmat*this->mpswf)<<endl;
    cout<<"<P2>="<<this->mpswf.dot(this->Invmat*this->Invmat*this->mpswf)<<endl;
    cout<<"<PH> ="<<this->mpswf.dot(this->PHmat*this->mpswf)<<endl;
    cout<<"<PH2>="<<this->mpswf.dot(this->PHmat*this->Invmat*this->mpswf)<<endl;
    this->print_mpswf();
    
    cout<<"mpswf1.mpswf2="<<mpswf1.dot(mpswf2)<<endl;
    
    cout<<"mpswf1.Invmat.mpswf2="<<mpswf1.dot(this->Invmat*mpswf2)<<endl;
    cout<<"mpswf1.PHmat.mpswf2 ="<<mpswf1.dot(this->PHmat*mpswf2)<<endl;
}
void Scars::print_mpswf(){
    cout<<"@ print mpswf @"<<endl;
    for (int i=0; i<this->statelist_M.size(); i++) {
        for (int j=0; j<this->No; j++) {
            cout<<setw(2)<<this->statelist_M[i].obasis[j]<<" ";
        }
        cout<<",   "<<setw(6)<<this->mpswf[i]<<endl;
    }
    cout<<"@ done  mpswf @"<<endl;
}
//void Scars::generatebasis_spin(int rangeS){
//    if (this->needM) {
//        this->generatebasis_b();
//    }
//    else {
//        this->generatebasis_b_noM();
//    }
//
//    //cout<<"in generatebasis_spin, needM="<<needM<<endl;
//    //this->print_statelist();
//    vector<state> statelist_M2; statelist_M2.clear();
//    bool keep=true;
//    for (int i=0; i<this->statelist_M.size(); i++) {
//        keep=true;
//        for (int j=0; j<this->statelist_M[i].obasis.size(); j++) {
//            //cout<<"i="<<i<<" j="<<j<<endl;
//            if (this->statelist_M[i].obasis[j]>rangeS) {
//                //cout<<"break"<<endl;
//                keep=false;
//                break;
//            }
//        }
//        if (keep) {
//            statelist_M2.push_back(this->statelist_M[i]);
//        }
//    }
//    this->statelist_M=statelist_M2;
//    //    this->print_statelist();
//
//    //    this->ini_hoplist();
//}
Scars::~Scars(){};

void quantum_scar(string filename){
    int no, rangeS, cutsite, nvec1, nvec2;
    string bound_cond;
    
    ifstream infile(filename);
    infile>>no>>rangeS;
    infile>>bound_cond;
    infile>>cutsite;
    infile>>nvec1>>nvec2;
    infile.close();
    
    cout<<"no, rangeS, boundary_cond="<<no<<" "<<rangeS<<" "<<bound_cond<<endl;
    
    Scars qscars(no, rangeS, bound_cond, nvec1, nvec2);
    
    cout<<"H.dim()="<<qscars.bitlist.size()<<endl;
    
    qscars.ee_setup(0, qscars.No/2, qscars.bitlist);
    
    bool calculatees=true, calculateph=true, showes=true, showph=true;
    qscars.diag_scar_H_inv("both", calculatees, calculateph);
    
    qscars.Show_Eigen_Sets(50001, showes, showph);
    qscars.Print_Eigen_Sets(50001, "E_S_Plot.txt", showes, showph);
    
//    cout<<"Hamiltonian matrix is\n"<<qscars.get_H()<<endl;
    
//    cout<<"\n\n Basis is"<<endl; qscars.print_statelist();
    
//    //qscars.print_statelist();
//    qscars.diag_scar();
//    qscars.show_scar_energy(20);
    
//    cout<<"Invmat=\n"<<qscars.Invmat<<endl;
//    cout<<"PHmat =\n"<<qscars.PHmat<<endl;
    
//    qscars.test_mps_wfs();
    
//    qscars.print_commutators();
    
//    qscars.print_statelist();
//    qscars.ee_setup(0, 3, qscars.bitlist);
//    qscars.print_trunc_states();
//
//    //TEST ENTANGLEMENT.
//    qscars.diag_scar();
//    qscars.show_Ground_State(true);
//    Eigen::MatrixXd rho2;
//    qscars.ee_compute_rho(Xcd_to_Xd(qscars.Ground_State), rho2, qscars.bitlist);
//    cout<<"entropy="<<qscars.ee_eval_rho(rho2)<<endl;
//
//    cout<<"rho2=\n"<<rho2<<endl;
}
