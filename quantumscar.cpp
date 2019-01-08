//
//  quantumscar.cpp
//  quantumscar
//
//  Created by Jie Wang on 12/10/18.
//  Copyright © 2018 Jie Wang. All rights reserved.
//

#include "quantumscar.h"

Scars::Scars(){};
Scars::Scars(int no_t, int ranges_t, string type_t, bool constrain_t, string bound_cond_t):rangeS(ranges_t), type(type_t), constrain(constrain_t), bound_cond(bound_cond_t) {
    this->init(no_t);
    this->out_info();
    
    this->generate_hilberspace();
    cout<<"done make hilbertspace"<<endl;
    
    this->setup_symmetry();
    cout<<"done setup symmetry"<<endl;
    
    this->diag_setup(this->type);
    
//    cout<<"$ start diagonalization $"<<endl;
//    //this->diag_scar();
//    this->diag_scar_H_inv("both");
////    this->diag_scar_H_spectra("both", this->nvec1, this->nvec2);
//    cout<<"$ done  diagonalization $"<<endl;
//
//    this->print_Eigen_Sets(200);
}
void Scars::out_info(){
    cout<<"no, rangeS, type, constrain, boundary_cond="<<this->No<<" "<<this->rangeS<<" "<<this->type<<" "<<this->constrain<<" "<<this->bound_cond<<endl;
}
void Scars::init(int no_t){
    this->statelist_M.clear(); this->matrixlist_M.clear(); this->one=(state_int) 1;
    //TODO: only for spin-half, it is "fermion".
    this->No=no_t; this->Np=0; this->momentum=0; this->geometry="sphere"; this->statistics="fermion"; this->shift=vector<double>(2, 0.); this->useqh=false; this->needM=false;
    
    if (this->No%2==0) {
        this->usefullsymmetry(true);
    }
    else {
        this->usefullsymmetry(false);
    }
}
void Scars::usefullsymmetry(bool use){
    if (use) {
        cout<<"@ Use Inversion 'P' and Translation symmetry 'T^(L/2)' @"<<endl;
        this->fullsymmetry=true;
    }
    else {
        cout<<"@ Use Inversion 'P' only @"<<endl;
        this->fullsymmetry=false;
    }
}
void Scars::setup_symmetry(){
    if (this->bound_cond=="PBC") {
        this->Translation_setup();
        cout<<"done translation"<<endl;
        this->Inversion_setup();
        cout<<"done inversion"<<endl;
        this->ParticleHole_setup();
        cout<<"done particlehole"<<endl;
    }
    else if (this->bound_cond=="OBC") {
        this->Inversion_setup();
        this->ParticleHole_setup();
    }
}
void Scars::diag_setup(string type_t){
    this->make_mlist(type_t);
    cout<<"done  make mlist"<<endl;
    
    this->H_spamatrix=Eigen::SparseMatrix<double>(this->bitlist.size(), this->bitlist.size());
    vector<Eigen::Triplet<double>> triplets=mlist_to_ftriplet(this->matrixlist_M);
    this->H_spamatrix.setFromTriplets(triplets.begin(), triplets.end());
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
void Scars::generate_hilberspace(){
    this->basis.clear();
    for (state_int i=0; i<this->one<<this->No; i++) {
        if (this->legal_state(i) and this->constrain) {
            this->basis.push_back(i);
        }
        else if (!this->constrain) {
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
void Scars::make_mlist(string type_t){
    this->matrixlist_M.clear();
    int N=this->bitlist.size();
    
    if (type_t=="PXP") {
        for (unsigned int i=0; i<N; i++) {
            matrix k;
            k.ket=i;
            state_int bit_tmp0=this->bitlist[i], bit_tmp;
            for (int j=0; j<this->No; j++) {
                
                bit_tmp=bit_tmp0 ^ this->one<<j;
                if (!this->legal_state(bit_tmp) and this->constrain) {
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
    else if (type_t=="PX2P") {
        for (unsigned i=0; i<N; i++) {
            matrix k;
            k.ket=i;
            state_int bit_tmp0=this->bitlist[i], bit_tmp;
            for (int j1=0; j1<this->No; j1++) {
                for (int j2=j1; j2<this->No; j2++) {
                    
                    bit_tmp=bit_tmp0 ^ this->one<<j1;
                    bit_tmp=bit_tmp  ^ this->one<<j2;
                    
                    if (!this->legal_state(bit_tmp) and this->constrain) {
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
    else if (type_t=="PXPXP") {
        this->matrixlist_M.clear();
        this->make_mlist("PXP");
        
        vector<Eigen::Triplet<double>> triplets=mlist_to_ftriplet(this->matrixlist_M);
        
        Eigen::SparseMatrix<double> spatmp=Eigen::SparseMatrix<double>(N, N), spatmp2;
        spatmp.setFromTriplets(triplets.begin(), triplets.end());
        spatmp2=spatmp*spatmp;
        
        this->matrixlist_M=triplet_to_mlist(to_triplets(spatmp2));
    }
    else if (type_t=="Pert") {
        double V; cout<<"type in 'V':"<<endl; cin>>V;
        cout<<"P.X.P + 1/V.(P.X.X.P- P.X.P.X.P);    V="<<V<<endl;
        
        vector<Eigen::Triplet<double>> triplets;
        Eigen::SparseMatrix<double> spa_PXP, spa_PX2P, spa_PXPXP, spa_Pert;
        
        spa_PXP=Eigen::SparseMatrix<double>(N, N); spa_PX2P=spa_PXP; spa_PXPXP=spa_PXP; spa_Pert=spa_PXP;
        
        
        this->make_mlist("PXP");
        triplets=mlist_to_ftriplet(this->matrixlist_M);
        spa_PXP.setFromTriplets(triplets.begin(), triplets.end());
        
        this->make_mlist("PX2P");
        triplets=mlist_to_ftriplet(this->matrixlist_M);
        spa_PX2P.setFromTriplets(triplets.begin(), triplets.end());
        
        this->make_mlist("PXPXP");
        triplets=mlist_to_ftriplet(this->matrixlist_M);
        spa_PXPXP.setFromTriplets(triplets.begin(), triplets.end());
        
        spa_Pert=spa_PXP+1./V*(spa_PX2P-spa_PXPXP);
        this->matrixlist_M=triplet_to_mlist(to_triplets(spa_Pert));
    }
    else {
        cout<<"unrecognized type"<<endl;
        exit(0);
    }
    
    this->matrixlist_M=mergemlist(sortmlist(this->matrixlist_M));
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
vector<double> Scars::get_energies(string type){
    vector<double> outputs;
    if (type=="pls") {
        for (unsigned i=0; i<this->Eigen_Sets_pls.size(); i++) {
            outputs.push_back(this->Eigen_Sets_pls[i].eval);
        }
    }
    else if (type=="min") {
        for (unsigned i=0; i<this->Eigen_Sets_min.size(); i++) {
            outputs.push_back(this->Eigen_Sets_min[i].eval);
        }
    }
    else if (type=="both") {
        for (unsigned i=0; i<this->Eigen_Sets.size(); i++) {
            outputs.push_back(this->Eigen_Sets[i].eval);
        }
    }
    else {
        cout<<"unrecognized type"<<endl;
        exit(0);
    }
    return outputs;
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
void Scars::make_Invmat(){
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
}
void Scars::Inversion_setup(){
    this->make_Invmat();
    this->make_inversion_proj();
}
void Scars::make_inversion_proj(){
    //setup the dimension of inversion symmetry, anti-symmetry sector.
    int dim_pls=0, dim_min=0;
    for (unsigned i=0; i<this->bitlist.size(); i++) {
        if (this->bitlist[i]==invert_bits(this->bitlist[i], this->No)) {
            dim_pls++;
        }
    }
    dim_pls=(this->bitlist.size()+dim_pls)/2;
    dim_min=this->bitlist.size()-dim_pls;
    this->P_pls=Eigen::MatrixXd::Zero(this->bitlist.size(), dim_pls);
    this->P_min=Eigen::MatrixXd::Zero(this->bitlist.size(), dim_min);
    
    //make inversion projectors;
    vector<state_int>::iterator binaryit;
    vector<int> found_states; found_states.clear();
    
    int counter_pls=0, counter_min=0;
    for (unsigned i=0; i<this->bitlist.size(); i++) {
        if (find(found_states.begin(), found_states.end(), i)!=found_states.end()) {
            continue;
        }
        state_int tmpbit=invert_bits(this->bitlist[i], this->No);
        if (this->bitlist[i]==tmpbit) {
            this->P_pls(i, counter_pls++)= 1.;
        }
        else {
            binaryit=lower_bound(this->bitlist.begin(), this->bitlist.end(), tmpbit);
            if (binaryit!=this->bitlist.end() and *binaryit==tmpbit) {
                int bitpos=binaryit-this->bitlist.begin();
                found_states.push_back(bitpos);
                this->P_pls(i,      counter_pls)   = 1./sqrt(2.);
                this->P_pls(bitpos, counter_pls++) = 1./sqrt(2.);
                this->P_min(i,      counter_min)   =-1./sqrt(2.);
                this->P_min(bitpos, counter_min++) = 1./sqrt(2.);
            }
            else {
                print_bit(this->bitlist[i], this->No); print_bit(tmpbit, this->No); cout<<"ind="<<binaryit-this->bitlist.begin()<<endl<<endl;
                cout<<"cannot find translated bit"<<endl;
                exit(0);
            }
        }
    }
}
void Scars::proj_H_invsector(){
    this->H_inv_pls=this->P_pls.adjoint()*this->H_spamatrix*this->P_pls;
    this->H_inv_min=this->P_min.adjoint()*this->H_spamatrix*this->P_min;
}
void Scars::diag_scar_H_inv(string invpls, bool calculatetrans, bool calculatees, bool calculateph){
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
                if (calculatetrans) {
                    Eigen::VectorXd vectmp=Xcd_to_Xd(this->P_pls*Xd_to_Xcd(eset_tmp.evec));
                    if (this->No%2==0) {
                        eset_tmp.trans=vectmp.dot(this->Transmat_halfL*vectmp);
                    }
                    else {
                        eset_tmp.trans=0.;
                    }
                }
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
        this->Eigen_Sets=this->Eigen_Sets_pls;
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
                if (calculatetrans) {
                    Eigen::VectorXd vectmp=Xcd_to_Xd(this->P_min*Xd_to_Xcd(eset_tmp.evec));
                    if (this->No%2==0) {
                        eset_tmp.trans=vectmp.dot(this->Transmat_halfL*vectmp);
                    }
                    else {
                        eset_tmp.trans=0;
                    }
                }
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
        this->Eigen_Sets=this->Eigen_Sets_min;
    }
    else if (invpls=="both") {
        this->diag_scar_H_inv("pls", calculatetrans, calculatees, calculateph);
        this->diag_scar_H_inv("min", calculatetrans, calculatees, calculateph);
        this->Eigen_Sets.clear();
        this->Eigen_Sets.insert(this->Eigen_Sets.end(), this->Eigen_Sets_pls.begin(), this->Eigen_Sets_pls.end());
        this->Eigen_Sets.insert(this->Eigen_Sets.end(), this->Eigen_Sets_min.begin(), this->Eigen_Sets_min.end());
    }
    
    sort(this->Eigen_Sets.begin(), this->Eigen_Sets.end(), compare_eigen_set);
}
bool Scars::diag_fullsymmetry(bool inv, bool trans, bool calculatees){
    if (trans) {
        if (!this->make_invtrans_proj(inv)) {
            cout<<"#false#"<<endl;
            return false;
        }
        else {
            //cout<<this->invtrans_proj.cols()<<endl;
            this->reducedH=this->invtrans_proj.adjoint()*this->H_spamatrix*this->invtrans_proj;
            
            //cout<<"@@ begin diagonalize (-) @@"<<endl;
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,-1,-1>> es;
            es.compute(this->reducedH);
            
            this->Eigen_Sets.clear();
            for (int i=0; i<es.eigenvalues().size(); i++) {
                eigen_set eset_tmp;
                eset_tmp.parity=inv;
                eset_tmp.trans=2*trans-1;
                eset_tmp.evec=es.eigenvectors().col(i);
                eset_tmp.eval=es.eigenvalues()[i];
                if (calculatees) {
                    Eigen::MatrixXd rho2;
                    this->ee_compute_rho(Xcd_to_Xd(this->invtrans_proj*Xd_to_Xcd(eset_tmp.evec)), rho2, this->bitlist);
                    eset_tmp.entanglement=this->ee_eval_rho(rho2);
                }
                this->Eigen_Sets.push_back(eset_tmp);
            }
        }
        sort(this->Eigen_Sets.begin(), this->Eigen_Sets.end(), compare_eigen_set);
        return true;
    }
    else {
        string invflag;
        if (inv) {
            invflag="pls";
        }
        else {
            invflag="min";
        }
        
        bool trans=false, calcualteph=false;
        this->diag_scar_H_inv(invflag, trans, calculatees, calcualteph);
        return true;
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
void Scars::Show_Eigen_Sets(int range, bool printtrans, bool printes, bool printph){
    range/=2;
    int N=this->Eigen_Sets.size();
    cout<<"N="<<N<<endl;
    
    cout<<"E,   Inv,   ";
    if (printtrans) {
        cout<<"Trans,   ";
    }
    if (printes) {
        cout<<"Entropy,   ";
    }
    if (printph) {
        cout<<"PH,   ";
    }
    cout<<endl;
    
    if (N%2==1) {
        N=(N-1)/2;
        range=min(range, N);
        for (int i=N-range; i<=N+range; i++) {
            cout<<setprecision(6)<<setw(15)<<chop(this->Eigen_Sets[i].eval)<<",   "<<printparity(this->Eigen_Sets[i].parity);
            if (printtrans) {
                cout<<",   "<<setw(6)<<chop(this->Eigen_Sets[i].trans);
            }
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
            if (printtrans) {
                cout<<",   "<<setw(6)<<chop(this->Eigen_Sets[i].trans);
            }
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
//    cout<<"Zero mode dimension = "<<zeromodedim<<endl;
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
void Scars::pick_Eigen_Sets(double E1, double E2, double S1, double S2, vector<eigen_set>& ret, vector<Eigen::VectorXd>& spectrum){
    ret.clear(); spectrum.clear();
    
    for (unsigned int i=0; i<this->Eigen_Sets.size(); i++) {
        if (this->Eigen_Sets[i].eval < E2 and this->Eigen_Sets[i].eval > E1 and this->Eigen_Sets[i].entanglement < S2 and this->Eigen_Sets[i].entanglement > S1) {
            ret.push_back(this->Eigen_Sets[i]);
        }
    }
    for (unsigned int i=0; i<ret.size(); i++) {
//        Eigen::MatrixXd rho2;
//        this->ee_compute_rho(Xcd_to_Xd(ret[i].evec), rho2, this->bitlist);
//        spectrum.push_back(this->espectrum_eval_rho(rho2));
        spectrum.push_back(this->get_espectrum(Xcd_to_Xd(ret[i].evec), this->bitlist));
    }
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
    if (this->No%2==0) {
        this->Transmat_halfL=Trans_pow(this->No/2);
    }
}
bool Scars::make_invtrans_proj(bool inv){
    Eigen::MatrixXd InvTransMat, invmatrix;
    if (inv) {
        invmatrix=this->P_pls;
        InvTransMat=this->P_pls.adjoint()*this->Transmat_halfL*this->P_pls;
    }
    else if (!inv) {
        invmatrix=this->P_min;
        InvTransMat=this->P_min.adjoint()*this->Transmat_halfL*this->P_min;
    }
    
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,-1,-1>> es;
    es.compute(InvTransMat);
    vector<double> InvTransEval=vector<double>(es.eigenvalues().size(), 0.);
    for (int i=0; i<InvTransEval.size(); i++) {
        InvTransEval[i]=es.eigenvalues()[i];
    }
    vector<Eigen::VectorXd> InvTransEvec=vector<Eigen::VectorXd>(es.eigenvalues().size());
    for (int i=0; i<InvTransEvec.size(); i++) {
        InvTransEvec[i]=Eigen::VectorXd::Zero(invmatrix.cols());
        for (int j=0; j<invmatrix.cols(); j++) {
            InvTransEvec[i](j)=es.eigenvectors().col(i)[j];
        }
    }
    int Eval_N=0;
    for (unsigned i=0; i<InvTransEval.size(); i++) {
        if (!(inv ^ this->one) and InvTransEval[i]>0) {
            Eval_N++;
        }
        else if ((inv ^ this->one) and InvTransEval[i]<0) {
            Eval_N++;
        }
    }
    
    this->invtrans_proj=Eigen::MatrixXd::Zero(this->bitlist.size(), Eval_N);
    int counter=0;
    for (unsigned i=0; i<InvTransEval.size(); i++) {
        if (inv and InvTransEval[i]>0) {
            this->invtrans_proj.col(counter++)=this->P_pls*InvTransEvec[i];
        }
        else if (!inv and InvTransEval[i]<0) {
            this->invtrans_proj.col(counter++)=this->P_min*InvTransEvec[i];
        }
    }
    //invtrans_proj.ajd() * invtrans_proj = Identity.
    
    if (invtrans_proj.cols()==0) {
        return false;
    }
    else {
        return true;
    }
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
        cout<<"[H, T]=\n"<<commutator(this->H_spamatrix, this->Transmat, true)<<endl;
        
        cout<<"[P, T]=\n"<<commutator(this->Invmat, this->Transmat, true)<<endl;
        
        if (this->No%2==0) {
            cout<<"[P, T^(L/2)]=\n"<<commutator(this->Invmat, this->Trans_pow(this->No/2), true)<<endl;
        }
        else {
            cout<<"[P, T^L]=\n"<<commutator(this->Invmat, this->Trans_pow(this->No), true)<<endl;
        }
        
        //cout<<"{P, T}=\n"<<commutator(this->Invmat, this->Transmat, false)<<endl;
        
        //cout<<"[P, T^L]=\n"<<commutator(this->Invmat, this->Trans_pow(this->No), true)<<endl;
        //cout<<"{P, T^L}=\n"<<commutator(this->Invmat, this->Trans_pow(this->No), false)<<endl;
        
        //cout<<"[P, T^(L/2)]=\n"<<commutator(this->Invmat, this->Trans_pow(this->No/2), true)<<endl;
        //cout<<"{P, T^(L/2)}=\n"<<commutator(this->Invmat, this->Trans_pow(this->No/2), false)<<endl;
        
        cout<<"[C, T]=\n"<<commutator(this->PHmat, this->Transmat, true)<<endl;
        cout<<"{C, T}=\n"<<commutator(this->PHmat, this->Transmat, false)<<endl;
        
        //cout<<"[C, T^L]=\n"<<commutator(this->PHmat, this->Trans_pow(this->No), true)<<endl;
        //cout<<"{C, T^L}=\n"<<commutator(this->PHmat, this->Trans_pow(this->No), false)<<endl;
        
        //cout<<"[C, T^(L/2)]=\n"<<commutator(this->PHmat, this->Trans_pow(this->No/2), true)<<endl;
        //cout<<"{C, T^(L/2)}=\n"<<commutator(this->PHmat, this->Trans_pow(this->No/2), false)<<endl;
        
        cout<<"T^L = \n"<<this->Trans_pow(this->No)<<endl;
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
        
        this->mpswf=Eigen::VectorXd::Zero(this->bitlist.size());
        
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
}
void Scars::mps_run(){
    if (this->No%2==1) {
        cout<<"L%2==1, quit"<<endl;
        exit(0);
    }
    
    if (this->bound_cond=="OBC") {
        B0=Eigen::MatrixXd(2, 3);
        B1=Eigen::MatrixXd(2, 3);
        C0=Eigen::MatrixXd(3, 2);
        C1=Eigen::MatrixXd(3, 2);
        Eigen::VectorXd v0=Eigen::VectorXd(2);
        Eigen::VectorXd v1=Eigen::VectorXd(2);
        v0<<1.,  1.;
        v1<<1., -1.;
        B0  <<  1,  0,  0,  0,       1,  0;
        B1  <<  0,  0,  0,  sqrt(2.),0,  sqrt(2.);
        C0  <<  0, -1,  1,  0,       0,  0;
        C1  <<sqrt(2.),0,0, 0,-sqrt(2.), 0;
        
        this->mpswf0=Eigen::VectorXd::Zero(this->bitlist.size());
        this->mpswf1=Eigen::VectorXd::Zero(this->bitlist.size());
        
        for (int i=0; i<this->bitlist.size(); i++) {
            Eigen::MatrixXd matrixprod0=v0.transpose();
            Eigen::MatrixXd matrixprod1=v1.transpose();
            for (int j=0; j<this->No/2; j++) {
                if (this->statelist_M[i].obasis[2*j]==0) {
                    matrixprod0*=B0;
                    matrixprod1*=B0;
                }
                else {
                    matrixprod0*=B1;
                    matrixprod1*=B1;
                }
                if (this->statelist_M[i].obasis[2*j+1]==0) {
                    matrixprod0*=C0;
                    matrixprod1*=C0;
                }
                else {
                    matrixprod0*=C1;
                    matrixprod1*=C1;
                }
            }
            matrixprod0*=v1;
            matrixprod1*=v0;
            if (matrixprod0.size()!=1 or matrixprod1.size()!=1) {
                cout<<"matrixprod.size()!=1"<<endl;
                exit(0);
            }
            this->mpswf0(i)=matrixprod0(0, 0);
            this->mpswf1(i)=matrixprod1(0, 0);
        }
    }
    else {
        cout<<"cannot do PBC so far"<<endl;
        exit(0);
    }
    
    this->mpswf0/=this->mpswf0.norm();
    this->mpswf1/=this->mpswf1.norm();
    
    double P00=this->mpswf0.dot(this->Invmat*this->mpswf0);
    double P01=this->mpswf0.dot(this->Invmat*this->mpswf1);
    double P10=this->mpswf1.dot(this->Invmat*this->mpswf0);
    double P11=this->mpswf1.dot(this->Invmat*this->mpswf1);
//    cout<<"P00, P01, P10, P11=\n"<<P00<<" "<<P01<<" "<<P10<<" "<<P11<<endl;
    Eigen::MatrixXd MPS_Pmat_01=Eigen::MatrixXd::Zero(2,2);
    MPS_Pmat_01<<P00, P01, P10, P11;
    cout<<"MPS_Pmat_01=\n"<<MPS_Pmat_01<<endl<<" MPS_Pmat*MPS_Pmat=\n"<<MPS_Pmat_01*MPS_Pmat_01<<endl;
    
    cout<<"mps0*mps1="<<mpswf0.dot(mpswf1)<<endl;
    
    cout<<"\n\nL="<<this->No<<endl;
    cout<<"\n\n@@@ Info of MPS0"<<endl;
    cout<<"MPS0.spectrum=\n"<<this->get_espectrum(this->mpswf0, this->bitlist)<<endl;
    cout<<"MPS0.entropy="<<this->get_entropy(this->mpswf0, this->bitlist)<<endl;
    cout<<"<P>="<<chop(this->mpswf0.dot(this->Invmat*this->mpswf0))<<endl;
    cout<<"<H>="<<chop(this->mpswf0.dot(this->H_spamatrix*this->mpswf0))<<endl;
    cout<<"<H^2>-<H>^2="<<chop(this->mpswf0.dot(this->H_spamatrix*this->H_spamatrix*this->mpswf0)-pow(chop(this->mpswf0.dot(this->H_spamatrix*this->mpswf0)),2))<<endl;
    
    cout<<"\n\n@@@ Info of MPS1"<<endl;
    cout<<"MPS1.spectrum=\n"<<this->get_espectrum(this->mpswf1, this->bitlist)<<endl;
    cout<<"MPS1.entropy="<<this->get_entropy(this->mpswf1, this->bitlist)<<endl;
    cout<<"<P>="<<chop(this->mpswf1.dot(this->Invmat*this->mpswf1))<<endl;
    cout<<"<H>="<<chop(this->mpswf1.dot(this->H_spamatrix*this->mpswf1))<<endl;
    cout<<"<H^2>-<H>^2="<<chop(this->mpswf1.dot(this->H_spamatrix*this->H_spamatrix*this->mpswf1)-pow(chop(this->mpswf1.dot(this->H_spamatrix*this->mpswf1)),2))<<endl;
    
//    cout<<"\n\n"<<endl;
//    this->mpswf_pls=this->mpswf0+this->mpswf1;
//    this->mpswf_min=this->mpswf0-this->mpswf1;
//    this->mpswf_pls/=this->mpswf_pls.norm();
//    this->mpswf_min/=this->mpswf_min.norm();
//    cout<<"mpspls*Invmat*mpspls=\n"<<this->mpswf_pls.dot(this->Invmat*this->mpswf_pls)<<endl;
//    double Ppp=this->mpswf_pls.dot(this->Invmat*this->mpswf_pls);
//    double Ppm=this->mpswf_pls.dot(this->Invmat*this->mpswf_min);
//    double Pmp=this->mpswf_min.dot(this->Invmat*this->mpswf_pls);
//    double Pmm=this->mpswf_min.dot(this->Invmat*this->mpswf_min);
//    cout<<"Ppp, Ppm, Pmp, Pmm=\n"<<Ppp<<" "<<Ppm<<" "<<Pmp<<" "<<Pmm<<endl;
//    Eigen::MatrixXd MPS_Pmat_pm=Eigen::MatrixXd::Zero(2,2);
//    MPS_Pmat_pm<<Ppp, Ppm, Pmp, Pmm;
//    cout<<"MPS_Pmat_pm=\n"<<MPS_Pmat_pm<<endl<<" MPS_Pmat_pm*MPS_Pmat_pm=\n"<<MPS_Pmat_pm*MPS_Pmat_pm<<endl;
//    cout<<"mpsp*mpsm="<<mpswf_pls.dot(mpswf_min)<<endl;
//
//    cout<<"MPS_pls.spectrum=\n"<<this->get_espectrum(this->mpswf_pls, this->bitlist)<<endl;
//    cout<<"MPS_pls.entropy=\n"<<this->get_entropy(this->mpswf_pls, this->bitlist)<<endl;
//
//    cout<<"MPS_pls*MPS_pls=\n"<<chop(this->mpswf_pls.dot(this->mpswf_pls))<<endl;
//    cout<<"MPS_pls*Inv*MPS_pls=\n"<<chop(this->mpswf_pls.dot(this->Invmat*this->mpswf_pls))<<endl;
//
//    cout<<"<H>=\n"<<chop(this->mpswf_pls.dot(this->H_spamatrix*this->mpswf_pls))<<endl;
//    cout<<"<H^2>-<H>^2=\n"<<chop(this->mpswf_pls.dot(this->H_spamatrix*this->H_spamatrix*this->mpswf_pls)-pow(chop(this->mpswf_pls.dot(this->H_spamatrix*this->mpswf_pls)), 2))<<endl;
}
void Scars::test_mps_wfs(){
    vector<Eigen::VectorXd> entanglement;
    
    this->mps_run(0);
    Eigen::VectorXd mpswf1=this->mpswf;
    cout<<"<P> ="<<this->mpswf.dot(this->Invmat*this->mpswf)<<endl;
    cout<<"<P2>="<<this->mpswf.dot(this->Invmat*this->Invmat*this->mpswf)<<endl;
    cout<<"<PH> ="<<this->mpswf.dot(this->PHmat*this->mpswf)<<endl;
    cout<<"<PH2>="<<this->mpswf.dot(this->PHmat*this->Invmat*this->mpswf)<<endl;
    cout<<"<H> ="<<chop(this->mpswf.transpose()*this->H_spamatrix*this->mpswf)<<endl;
    cout<<"<H2>="<<chop(this->mpswf.transpose()*this->H_spamatrix*this->H_spamatrix*this->mpswf)<<endl;
    entanglement=this->mps_entanglement();
    for (int i=0; i<entanglement.size(); i++) {
        double ent=0.;
        for (int j=0; j<entanglement[i].size(); j++) {
            if (abs(entanglement[i](j)>0)) {
                ent+=entanglement[i](j)*log(entanglement[i](j));
            }
        }
        cout<<"entropy="<<-ent<<endl;
    }
    //this->print_mpswf();
    
    this->mps_run(1);
    Eigen::VectorXd mpswf2=this->mpswf;
    cout<<"<P> ="<<this->mpswf.dot(this->Invmat*this->mpswf)<<endl;
    cout<<"<P2>="<<this->mpswf.dot(this->Invmat*this->Invmat*this->mpswf)<<endl;
    cout<<"<PH> ="<<this->mpswf.dot(this->PHmat*this->mpswf)<<endl;
    cout<<"<PH2>="<<this->mpswf.dot(this->PHmat*this->Invmat*this->mpswf)<<endl;
    cout<<"<H> ="<<chop(this->mpswf.transpose()*this->H_spamatrix*this->mpswf)<<endl;
    cout<<"<H2>="<<chop(this->mpswf.transpose()*this->H_spamatrix*this->H_spamatrix*this->mpswf)<<endl;
    entanglement=this->mps_entanglement();
    for (int i=0; i<entanglement.size(); i++) {
        double ent=0.;
        for (int j=0; j<entanglement[i].size(); j++) {
            if (abs(entanglement[i](j)>0)) {
                ent+=entanglement[i](j)*log(entanglement[i](j));
            }
        }
        cout<<"entropy="<<-ent<<endl;
    }
    //this->print_mpswf();
    
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
vector<Eigen::VectorXd> Scars::mps_entanglement(){
    vector<Eigen::VectorXd> entanglement; entanglement.clear();
    for (int i=0; i<2; i++) {
        cout<<"@ MPS STATE ["<<i<<"] @"<<endl;
        this->mps_run(i);
        Eigen::MatrixXd rho2;
        this->ee_compute_rho(Xcd_to_Xd(this->mpswf), rho2, this->bitlist);
        entanglement.push_back(this->espectrum_eval_rho(rho2));
    }
    return entanglement;
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
    int no, rangeS;
    bool constrain;
    string type, bound_cond;
    
    ifstream infile(filename);
    infile>>no>>rangeS;
    infile>>type>>constrain>>bound_cond;
    infile.close();
    
//    cout<<"no, rangeS, type, boundary_cond="<<no<<" "<<rangeS<<" "<<type<<" "<<bound_cond<<endl;
    
    Scars qscars(no, rangeS, type, constrain, bound_cond);
    
    cout<<"$ project H into inverse sector"<<endl;
    qscars.proj_H_invsector();
    
    cout<<"H.dim()="<<qscars.bitlist.size()<<endl;
    
//    qscars.print_commutators();
    
    qscars.ee_setup(0, qscars.No/2, qscars.bitlist);
    
    bool calculatetrans=true, calculatees=true, calculateph=true, showes=true, showph=false, showtrans=true;
//    qscars.diag_scar_H_inv("both", calculatetrans, calculatees, calculateph);
//    qscars.Show_Eigen_Sets(50000, showtrans, showes, showph);
    
    cout<<"\nFull Symmetry, true, true"<<endl;
    bool inv, usetrans;
    if (qscars.No%4==2) {
        inv=true;
    }
    else if (qscars.No%4==0) {
        inv=false;
    }
    else {
        cout<<"No odd, quit"<<endl;
        exit(0);
    }
    if (bound_cond=="PBC") {
        usetrans=true;
    }
    else {
        usetrans=false;
    }
    if (qscars.diag_fullsymmetry(inv, usetrans, calculatees)) {
        qscars.Show_Eigen_Sets(50000, showtrans, showes, showph);
        qscars.Print_Eigen_Sets(50000, "OBC/E_S_Plot_16.txt", showes, showph);
    }
    
    
//    qscars.diag_scar_H_inv("both", false, calculatees, calculateph);
//    qscars.Show_Eigen_Sets(50000, showtrans, showes, showph);
//
//    qscars.mps_run();
//
//
//    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,-1,-1>> es;
//    es.compute(qscars.get_H());
//    for (int i=0; i<es.eigenvalues().size(); i++) {
//        cout<<es.eigenvalues()(i)<<endl;
//    }
    
    
    
    vector<eigen_set> eigensets; vector<Eigen::VectorXd> spectrum;
    qscars.pick_Eigen_Sets(1.0, 5.0 ,0.4 ,1.0 , eigensets, spectrum);
    if (spectrum.size()!=1) {
        cout<<"spectrum.size()!=1, ="<<spectrum.size()<<endl;
//        exit(0);
    }
    for (int i=0; i<eigensets.size(); i++) {
        cout<<"i="<<i<<" eng="<<eigensets[i].eval<<" ent="<<eigensets[i].entanglement<<endl;
        
        vector<double> states; Xd_to_Vec(eigensets[i].evec, states);
        vector<double> spec2=qscars.entanglement_spectrum_SVD(states, qscars.bitlist, 1, 0);
        cout<<"spec2.size()="<<spec2.size()<<endl;
        for (int i=0; i<spec2.size(); i++) {
            cout<<spec2[i]<<endl;
        }
    }
    for (int i=0; i<spectrum.size(); i++) {
        
        
//        Eigen::VectorXd spectrums=this->get_espectrum(this->mpswf1, this->bitlist);
        
        
        
        ofstream outfilees("OBC/scar_spectrum_"+to_string((long long int)(no))+"_"+to_string((long long int)(i))+".txt");
        for (int j=0; j<spectrum[i].size(); j++) {
            outfilees<<spectrum[i](j)<<endl;
        }
        outfilees.close();
    }
    
    qscars.mps_run();
    
    
    
//    qscars.test_mps_wfs();
    
//    cout<<"\nFull Symmetry, true, false"<<endl;
//    if (qscars.diag_fullsymmetry(true, false, true)) {
//        qscars.Show_Eigen_Sets(50000, showtrans, showes, showph);
//    }
//
//    cout<<"\nFull Symmetry, false, true"<<endl;
//    if (qscars.diag_fullsymmetry(false, true, true)) {
//        qscars.Show_Eigen_Sets(50000, showtrans, showes, showph);
//    }
//
//    cout<<"\nFull Symmetry, false, false"<<endl;
//    if (qscars.diag_fullsymmetry(false, false, true)) {
//        qscars.Show_Eigen_Sets(50000, showtrans, showes, showph);
//    }
    
    
    
//    qscars.Print_Eigen_Sets(50000, "E_S_Plot.txt", showes, showph);

//    vector<double> leveldiff=leveldifference(qscars.get_energies("both"));
//
//    ofstream outfile("levelsta.txt");
//    for (int i=0; i<leveldiff.size(); i++) {
//        outfile<<leveldiff[i]<<endl;
//    }
//    outfile.close();
    
    
    
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
void quantum_scar_mps(string filename){
    int no, rangeS;
    bool constrain;
    string type, bound_cond;
    
    ifstream infile(filename);
    infile>>no>>rangeS;
    infile>>type>>constrain>>bound_cond;
    infile.close();
    
    Scars qscars(no, rangeS, type, constrain, bound_cond);
    
    qscars.ee_setup(0, qscars.No/2, qscars.bitlist);
    vector<Eigen::VectorXd> entanglement=qscars.mps_entanglement();
    for (int i=0; i<entanglement.size(); i++) {
        cout<<"i="<<endl;
        cout<<entanglement[i]<<endl;
        cout<<endl;
        
        double ent=0.;
        for (int j=0; j<entanglement[i].size(); j++) {
            if (abs(entanglement[i](j)>0)) {
                ent+=entanglement[i](j)*log(entanglement[i](j));
            }
        }
        cout<<"entropy="<<-ent<<endl;
    }
}
