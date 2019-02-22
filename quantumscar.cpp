//
//  quantumscar.cpp
//  quantumscar
//
//  Created by Jie Wang on 12/10/18.
//  Copyright © 2018 Jie Wang. All rights reserved.
//

#include "quantumscar.h"

Eigen::SparseMatrix<double> Scars::SparseRemoveRow(Eigen::SparseMatrix<double> input, int ind, vector<state_int>& bit_list, vector<state>& state_list) {
    //cout<<"input.rows()="<<input.rows()<<endl;
    
    if (ind<input.rows()) {
        bit_list.erase(bit_list.begin()+ind);
        state_list.erase(state_list.begin()+ind);
        
        Eigen::SparseMatrix<double> temp=Eigen::SparseMatrix<double>(input.cols(), input.rows()-1);
        
        for (int i=0; i<input.rows(); i++) {
            if (i<ind) {
                temp.col(i)=input.transpose().col(i);
            }
            else if (i>ind) {
                temp.col(i-1)=input.transpose().col(i);
            }
        }
        return temp.transpose();
    }
    else {
        cout<<"Cannot remove row."<<endl;
        exit(0);
    }
}
Eigen::SparseMatrix<double> SparseLinearComb(Eigen::SparseMatrix<double> input, int ind1, int ind2, double coeff1, double coeff2) {
    Eigen::SparseMatrix<double> temp=input.transpose();
    //temp.col(ind1)=temp.col(ind1)*coeff1+this->Invmat.transpose()*temp.col(ind2)*coeff2;
    temp.col(ind1)=temp.col(ind1)*coeff1+temp.col(ind2)*coeff2;
    if (temp.col(ind1).norm()!=0.) {
        temp.col(ind1)/=temp.col(ind1).norm();
    }
    //temp=SparseRemoveRow(temp.transpose(), ind2, bit_list, state_list);
    return temp.transpose();
}

Scars::Scars(){};
Scars::Scars(int no_t, int ranges_t, string type_t, bool constrain_t, string bound_cond_t, int replu_range_t, int connectX_t, bool quickmode_t):rangeS(ranges_t), type(type_t), constrain(constrain_t), bound_cond(bound_cond_t), replu_range(replu_range_t), connectX(connectX_t), quickmode(quickmode_t) {
    this->init(no_t);
    //this->out_info();
    
    this->generate_hilberspace();
    //cout<<"done make hilbertspace"<<endl;
    
    //this->print_statelist();
    
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
    cout<<"replusion_range, connected_X = "<<this->replu_range<<" "<<this->connectX<<endl;
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
    
    cout<<"L="<<this->No<<", type="<<this->type<<", "<<this->replu_range<<" "<<connectX<<endl;
}
void Scars::usefullsymmetry(bool use){
    if (use) {
        //cout<<"@ Use Inversion 'P' and Translation symmetry 'T^(L/2)' @"<<endl;
        this->fullsymmetry=true;
    }
    else {
        //cout<<"@ Use Inversion 'P' only @"<<endl;
        this->fullsymmetry=false;
    }
}
void Scars::setup_symmetry(){
    if (this->bound_cond=="PBC" and this->quickmode) {
        this->Translation_setup();
        //cout<<"done translation"<<endl;
    }
    else if (this->bound_cond=="PBC" and !this->quickmode) {
        this->Translation_setup();
        //cout<<"done translation"<<endl;
        this->Inversion_setup();
        //cout<<"done inversion"<<endl;
        this->make_inversion_proj();
        //cout<<"done inversion proj"<<endl;
        this->ParticleHole_setup();
        //cout<<"done particlehole"<<endl;
    }
    else if (this->bound_cond=="OBC") {
        this->Inversion_setup();
        this->ParticleHole_setup();
    }
}
void Scars::diag_setup(string type_t){
    this->make_mlist(type_t);
    //cout<<"done  make mlist"<<endl;
    
    this->H_spamatrix=Eigen::SparseMatrix<double>(this->bitlist.size(), this->bitlist.size());
    //cout<<"mlist_to_ftriplets"<<endl;
    vector<Eigen::Triplet<double>> triplets=mlist_to_ftriplet(this->matrixlist_M);
    //cout<<"set from triplets"<<endl;
    this->H_spamatrix.setFromTriplets(triplets.begin(), triplets.end());
    //cout<<"done diag setup"<<endl;
}
//inline bool Scars::legal_state(state_int input){
//    if (bound_cond=="PBC") {
//        for (int i=0; i<this->No; i++) {
//            if (this->one<<supermod(i, this->No)&input and this->one<<supermod(i+1, this->No)&input) {
//                return false;
//            }
//        }
//        return true;
//    }
//    else if (bound_cond=="OBC") {
//        //(i,i+i) w/ i<=No-2.
//        for (int i=0; i<this->No-1; i++) {
//            if (this->one<<supermod(i, this->No)&input and this->one<<supermod(i+1, this->No)&input) {
//                return false;
//            }
//        }
//        return true;
//    }
//    else {
//        cout<<"unrecognized boundary condition"<<endl;
//        exit(0);
//    }
//}
inline bool Scars::legal_state(state_int input){
    if (bound_cond=="PBC") {
        for (int i=0; i<this->No; i++) {
            if (this->one<<supermod(i, this->No)&input) {
                for (int j=1; j<=this->replu_range; j++) {
                    if (this->one<<supermod(i+j, this->No)&input) {
                        return false;
                    }
                }
            }
        }
        return true;
    }
    else if (bound_cond=="OBC") {
        //(i,i+i) w/ i<=No-2.
        for (int i=0; i<this->No; i++) {
            if (this->one<<supermod(i, this->No)&input) {
                for (int j=1; j<=this->replu_range; j++) {
                    if (this->one<<min(i+j, this->No)&input) {
                        return false;
                    }
//                    if (this->one<<supermod(i+j, this->No)&input) {
//                        return false;
//                    }
                }
            }
        }
        return true;
    }
    else {
        cout<<"unrecognized boundary condition"<<endl;
        exit(0);
    }
}
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
    
    cout<<"done generate hilbert space, dim()="<<this->bitlist.size()<<endl;
}
void Scars::make_mlist(string type_t){
    this->matrixlist_M.clear();
    state_int N=this->bitlist.size();
    
    if (type_t=="PXP" and this->connectX==1) {
        for (state_int i=0; i<N; i++) {
            matrix k;
            k.ket=i;
            state_int bit_tmp0=this->bitlist[i], bit_tmp;
            for (state_int j=0; j<this->No; j++) {
                
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
    else if (type_t=="PXP" and this->connectX!=1) {
        for (unsigned int i=0; i<N; i++) {
            matrix k;
            k.ket=i;
            state_int bit_tmp0=this->bitlist[i], bit_tmp;
            if (this->bound_cond=="PBC") {
                for (int j=0; j<this->No; j++) {
                    
                    bit_tmp=bit_tmp0 ^ this->one<<j;
                    for (int ind=1; ind<this->connectX; ind++) {
                        bit_tmp=bit_tmp ^ this->one<<supermod(j+ind, this->No);
                    }
                    
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
            else if (this->bound_cond=="OBC") {
                for (int j=0; j<this->No-this->connectX+1; j++) {
                    
                    bit_tmp=bit_tmp0 ^ this->one<<j;
                    for (int ind=1; ind<this->connectX; ind++) {
                        bit_tmp=bit_tmp ^ this->one<<j+ind;
                    }
                    
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
    }
    else if (type_t=="PX2P") {
        for (unsigned int i=0; i<N; i++) {
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
void Scars::diag_scar_H(){
    cout<<"@@ begin diagonalize @@, size="<<this->H_spamatrix.cols()<<endl;
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
    this->All_States=vector<Eigen::VectorXcd>(es.eigenvectors().cols(), Eigen::VectorXcd::Zero(this->Ground_State.size()));
    for (unsigned i=0; i<es.eigenvectors().cols(); i++) {
        for (unsigned j=0; j<this->Ground_State.size(); j++) {
            this->All_States[i](j)=es.eigenvectors().col(i)[j];
        }
    }
    cout<<"@@ done  diagonalize @@"<<endl;
}
int Scars::get_H_dim(){
    return this->bitlist.size();
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
void Scars::Inversion_setup(){
    this->Invmat=Eigen::SparseMatrix<double>(this->bitlist.size(), this->bitlist.size());
    
    vector<Eigen::Triplet<double>> triples; triples.clear();
    
    for (state_int i=0; i<this->bitlist.size(); i++) {
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
void Scars::make_inversion_proj(){
    //setup the dimension of inversion symmetry, anti-symmetry sector.
    int dim_pls=0, dim_min=0;
    for (state_int i=0; i<this->bitlist.size(); i++) {
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
    for (state_int i=0; i<this->bitlist.size(); i++) {
        if (i%1000==0) {
            cout<<"i="<<i<<endl;
        }
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
            
            int nthreads=omp_get_max_threads();
            vector<vector<eigen_set>> EigenSetspls(nthreads);
            for (int i=0; i<nthreads; i++) {
                EigenSetspls[i].clear();
            }
            this->Eigen_Sets_pls.clear();
            
            #pragma omp parallel for
            for (state_int i=0; i<es.eigenvalues().size(); i++) {
                state_int coren=omp_get_thread_num();
                eigen_set eset_tmp;
                eset_tmp.parity=true;
                eset_tmp.evec=this->P_pls*es.eigenvectors().col(i);
                eset_tmp.eval=es.eigenvalues()[i];
                eset_tmp.singletrans1=chop(eset_tmp.evec.dot(this->Transmat*eset_tmp.evec));
                eset_tmp.singletrans2=chop(eset_tmp.evec.dot(this->Transmat*this->Transmat*eset_tmp.evec));
                if (calculatetrans) {
                    if (this->No%2==0) {
                        eset_tmp.trans=eset_tmp.evec.dot(this->Transmat_halfL*eset_tmp.evec);
                    }
                    else {
                        eset_tmp.trans=0.;
                    }
                }
                if (calculatees) {
                    Eigen::MatrixXd rho2;
                    this->ee_compute_rho(eset_tmp.evec, rho2, this->bitlist);
                    eset_tmp.entanglement=this->ee_eval_rho(rho2);
                }
                if (calculateph) {
                    eset_tmp.particlehole=eset_tmp.evec.dot(this->PHmat*eset_tmp.evec);
                }
                EigenSetspls[coren].push_back(eset_tmp);
                //this->Eigen_Sets_pls.push_back(eset_tmp);
            }
            for (state_int i=0; i<nthreads; i++) {
                for (state_int j=0; j<EigenSetspls[i].size(); j++) {
                    this->Eigen_Sets_pls.push_back(EigenSetspls[i][j]);
                }
            }
            sort(this->Eigen_Sets_pls.begin(), this->Eigen_Sets_pls.end(), compare_eigen_set);
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
            
            int nthreads=omp_get_max_threads();
            vector<vector<eigen_set>> EigenSetsmin(nthreads);
            for (int i=0; i<nthreads; i++) {
                EigenSetsmin[i].clear();
            }
            this->Eigen_Sets_min.clear();
            
            #pragma omp parallel for
            for (int i=0; i<es.eigenvalues().size(); i++) {
                int coren=omp_get_thread_num();
                eigen_set eset_tmp;
                eset_tmp.parity=false;
                eset_tmp.evec=this->P_min*es.eigenvectors().col(i);
                eset_tmp.eval=es.eigenvalues()[i];
                eset_tmp.singletrans1=chop(eset_tmp.evec.dot(this->Transmat*eset_tmp.evec));
                eset_tmp.singletrans2=chop(eset_tmp.evec.dot(this->Transmat*this->Transmat*eset_tmp.evec));
                if (calculatetrans) {
                    if (this->No%2==0) {
                        eset_tmp.trans=eset_tmp.evec.dot(this->Transmat_halfL*eset_tmp.evec);
                    }
                    else {
                        eset_tmp.trans=0;
                    }
                }
                if (calculatees) {
                    Eigen::MatrixXd rho2;
                    this->ee_compute_rho(eset_tmp.evec, rho2, this->bitlist);
                    eset_tmp.entanglement=this->ee_eval_rho(rho2);
                }
                if (calculateph) {
                    eset_tmp.particlehole=eset_tmp.evec.dot(this->PHmat*eset_tmp.evec);
                }
                EigenSetsmin[coren].push_back(eset_tmp);
                //this->Eigen_Sets_min.push_back(eset_tmp);
            }
            for (int i=0; i<nthreads; i++) {
                for (int j=0; j<EigenSetsmin[i].size(); j++) {
                    this->Eigen_Sets_min.push_back(EigenSetsmin[i][j]);
                }
            }
            sort(this->Eigen_Sets_min.begin(), this->Eigen_Sets_min.end(), compare_eigen_set);
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
bool Scars::diag_scar_H_inv_halft(bool inv, bool trans, bool calculatees){
    bool zerosector=false;
    if (!this->make_invtrans_proj(inv, trans, zerosector)) {
        cout<<"#false#"<<endl;
        return false;
    }
    else {
        //cout<<this->invtrans_proj.cols()<<endl;
        this->reducedH=this->invtrans_proj.adjoint()*this->H_spamatrix*this->invtrans_proj;
        cout<<"diagfullsym,reducedh.dim="<<this->reducedH.cols()<<endl;
        
        //cout<<"@@ begin diagonalize (-) @@"<<endl;
        //cout<<"@@ begin diagonalize, reducedH.size="<<this->reducedH.cols()<<endl;
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,-1,-1>> es;
        es.compute(this->reducedH);
        
        int nthreads=omp_get_max_threads();
        vector<vector<eigen_set>> EigenSets(nthreads);
        for (int i=0; i<nthreads; i++) {
            EigenSets[i].clear();
        }
        this->Eigen_Sets.clear();
        
        #pragma omp parallel for
        for (int i=0; i<es.eigenvalues().size(); i++) {
            int coren=omp_get_thread_num();
            eigen_set eset_tmp;
            eset_tmp.parity=inv;
            eset_tmp.trans=2*trans-1;
            eset_tmp.evec=this->invtrans_proj*es.eigenvectors().col(i);
            eset_tmp.singletrans1=chop(eset_tmp.evec.dot(this->Transmat*eset_tmp.evec));
            eset_tmp.singletrans2=chop(eset_tmp.evec.dot(this->Transmat*this->Transmat*eset_tmp.evec));
            eset_tmp.eval=es.eigenvalues()[i];
            if (calculatees) {
                Eigen::MatrixXd rho2;
                this->ee_compute_rho(eset_tmp.evec, rho2, this->bitlist);
                eset_tmp.entanglement=this->ee_eval_rho(rho2);
            }
            EigenSets[coren].push_back(eset_tmp);
            //this->Eigen_Sets.push_back(eset_tmp);
        }
        for (int i=0; i<nthreads; i++) {
            for (int j=0; j<EigenSets[i].size(); j++) {
                this->Eigen_Sets.push_back(EigenSets[i][j]);
            }
        }
    }
    sort(this->Eigen_Sets.begin(), this->Eigen_Sets.end(), compare_eigen_set);
    return true;
//    if (trans) {
//        bool zerosector=false;
//        if (!this->make_invtrans_proj(inv, true, zerosector)) {
//            cout<<"#false#"<<endl;
//            return false;
//        }
//        else {
//            //cout<<this->invtrans_proj.cols()<<endl;
//            this->reducedH=this->invtrans_proj.adjoint()*this->H_spamatrix*this->invtrans_proj;
//            cout<<"diagfullsym,reducedh.dim="<<this->reducedH.cols()<<endl;
//
//            //cout<<"@@ begin diagonalize (-) @@"<<endl;
//            //cout<<"@@ begin diagonalize, reducedH.size="<<this->reducedH.cols()<<endl;
//            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,-1,-1>> es;
//            es.compute(this->reducedH);
//
//            int nthreads=omp_get_max_threads();
//            vector<vector<eigen_set>> EigenSets(nthreads);
//            for (int i=0; i<nthreads; i++) {
//                EigenSets[i].clear();
//            }
//            this->Eigen_Sets.clear();
//
//            #pragma omp parallel for
//            for (int i=0; i<es.eigenvalues().size(); i++) {
//                int coren=omp_get_thread_num();
//                eigen_set eset_tmp;
//                eset_tmp.parity=inv;
//                eset_tmp.trans=2*trans-1;
//                eset_tmp.evec=this->invtrans_proj*es.eigenvectors().col(i);
//                eset_tmp.singletrans1=chop(eset_tmp.evec.dot(this->Transmat*eset_tmp.evec));
//                eset_tmp.singletrans2=chop(eset_tmp.evec.dot(this->Transmat*this->Transmat*eset_tmp.evec));
//                eset_tmp.eval=es.eigenvalues()[i];
//                if (calculatees) {
//                    Eigen::MatrixXd rho2;
//                    this->ee_compute_rho(eset_tmp.evec, rho2, this->bitlist);
//                    eset_tmp.entanglement=this->ee_eval_rho(rho2);
//                }
//                EigenSets[coren].push_back(eset_tmp);
//                //this->Eigen_Sets.push_back(eset_tmp);
//            }
//            for (int i=0; i<nthreads; i++) {
//                for (int j=0; j<EigenSets[i].size(); j++) {
//                    this->Eigen_Sets.push_back(EigenSets[i][j]);
//                }
//            }
//        }
//        sort(this->Eigen_Sets.begin(), this->Eigen_Sets.end(), compare_eigen_set);
//        return true;
//    }
//    else {
//        string invflag;
//        if (inv) {
//            invflag="pls";
//        }
//        else {
//            invflag="min";
//        }
//
//        bool trans=false, calcualteph=false;
//        this->diag_scar_H_inv(invflag, trans, calculatees, calcualteph);
//        return true;
//    }
}
bool Scars::diag_scar(bool inv, int Ky, string mode, bool calculatees, bool calculatet, bool calculatett, bool calculatethalf) {
    bool trans;
    if (Ky==1) {
        trans=true;
    }
    else if (Ky==-1) {
        trans=false;
    }
    else {
        cout<<"wrong Ky"<<endl;
        exit(0);
    }
    
    Eigen::SparseMatrix<double> localshrink;
    if (mode=="t") {
        this->makeScarShrinker_trans(Ky);
        this->makeScarShrinker_trans_inv(inv);
        localshrink=this->shrinkMatrix_trans_inv;
    }
    else if (mode=="halft") {
        this->make_invtrans_proj(inv, trans);
        localshrink=this->invtrans_proj.adjoint();
    }
    
    if (localshrink.rows()==0) {
        return false;
    }
    
    this->reducedH=localshrink*this->H_spamatrix*localshrink.adjoint();
    cout<<"reducedH.dim()="<<this->reducedH.cols()<<endl;
    
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,-1,-1>> es;
    es.compute(this->reducedH);
    cout<<"done diagonalization"<<endl;
    
    int nthreads=omp_get_max_threads();
    vector<vector<eigen_set>> EigenSets(nthreads);
    for (int i=0; i<nthreads; i++) {
        EigenSets[i].clear();
    }
    this->Eigen_Sets.clear();
    
    #pragma omp parallel for
    for (state_int i=0; i<es.eigenvalues().size(); i++) {
        state_int coren=omp_get_thread_num();
        eigen_set eset_tmp;
        eset_tmp.parity=inv;
        //eset_tmp.trans=2*trans-1;
        eset_tmp.evec=localshrink.adjoint()*es.eigenvectors().col(i);
        eset_tmp.eval=es.eigenvalues()[i];
        
        //cout<<"E="<<eset_tmp.eval<<endl;
        if (calculatethalf) {
            eset_tmp.trans=chop(eset_tmp.evec.dot(this->Transmat_halfL*eset_tmp.evec));
        }
        if (calculatet) {
            eset_tmp.singletrans1=chop(eset_tmp.evec.dot(this->Transmat*eset_tmp.evec));
        }
        if (calculatett) {
            eset_tmp.singletrans2=chop(eset_tmp.evec.dot(this->Transmat*this->Transmat*eset_tmp.evec));
        }
        if (calculatees) {
            Eigen::MatrixXd rho2;
            Eigen::MatrixXd svdmat;
            //cout<<"to redmindmat"<<endl;
            //this->ee_compute_rho_redmindmat(eset_tmp.evec, rho2);
            //cout<<"done redmindmat"<<endl;
            //cout<<"to compute ent, dm"<<endl;
            //eset_tmp.entanglement=this->ee_eval_rho(rho2);
            //cout<<"done ent, dm, "<<eset_tmp.entanglement<<endl;
            
            //cout<<"to svd"<<endl;
            this->generate_svdmat(eset_tmp.evec, svdmat);
            //cout<<"done svd"<<endl;
            //cout<<"to compute end, svd"<<endl;
            eset_tmp.entanglement=ee_eval_bdcsvd(svdmat);
            //cout<<"done ent, svd, "<<eset_tmp.entanglement<<endl;
            //cout<<"done calculate ent"<<endl;
        }
        
        EigenSets[coren].push_back(eset_tmp);
        //this->Eigen_Sets.push_back(eset_tmp);
    }
    for (state_int i=0; i<nthreads; i++) {
        for (state_int j=0; j<EigenSets[i].size(); j++) {
            this->Eigen_Sets.push_back(EigenSets[i][j]);
        }
    }
    
    //cout<<"shrinkermatrix=\n"<<localshrink<<endl;
    //cout<<"to sort"<<endl;
    //sort(this->Eigen_Sets.begin(), this->Eigen_Sets.end(), compare_eigen_set);
    return true;
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
    //cout<<"N="<<N<<endl;
    
    //cout<<"statelist=\n"; this->print_statelist();
    
    cout<<"      E,   Inv,   ";
    if (printtrans) {
        cout<<"t^{L/2},   ";
        cout<<"t,   ";
        cout<<"t*t,   ";
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
            cout<<setprecision(6)<<setw(10)<<chop(this->Eigen_Sets[i].eval)<<",   "<<printparity(this->Eigen_Sets[i].parity);
            if (printtrans) {
                cout<<",   "<<setw(6)<<chop(this->Eigen_Sets[i].trans);
                cout<<",   "<<setprecision(5)<<setw(10)<<chop(this->Eigen_Sets[i].singletrans1);
                cout<<",   "<<setprecision(5)<<setw(10)<<chop(this->Eigen_Sets[i].singletrans2);
            }
            if (printph) {
                cout<<",   "<<setw(6)<<chop(this->Eigen_Sets[i].particlehole);
            }
            if (printes) {
                cout<<",   "<<setw(15)<<chop(this->Eigen_Sets[i].entanglement);
            }
            //cout<<endl<<this->Eigen_Sets[i].evec<<endl;
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
                cout<<",   "<<setprecision(5)<<setw(10)<<chop(this->Eigen_Sets[i].singletrans1);
                cout<<",   "<<setprecision(5)<<setw(10)<<chop(this->Eigen_Sets[i].singletrans2);
            }
            if (printph) {
                cout<<",   "<<setw(6)<<chop(this->Eigen_Sets[i].particlehole);
            }
            if (printes) {
                cout<<",   "<<setw(15)<<chop(this->Eigen_Sets[i].entanglement);
            }
            //cout<<endl<<this->Eigen_Sets[i].evec<<endl;
            cout<<endl;
        }
    }
    
    int zeromodedim=0;
    for (int i=0; i<this->Eigen_Sets.size(); i++) {
        if (abs(this->Eigen_Sets[i].eval)<1e-12) {
            zeromodedim++;
        }
    }
    //cout<<"Zero mode dimension = "<<zeromodedim<<endl;
}
void Scars::Print_Eigen_Sets2(int range, string filename, bool printes, bool printtrans){
    ofstream outfile("data/"+filename);
    
    outfile<<"      E,      Inv,      ";
    if (printtrans) {
        outfile<<"t^{L/2},      ";
        outfile<<"t,        ";
        outfile<<"t*t,         ";
    }
    if (printes) {
        outfile<<"Entropy,      ";
    }
    outfile<<endl;
    
    if (range!=-1) {
        range/=2;
        int N=this->Eigen_Sets.size();
        if (N%2==1) {
            N=(N-1)/2;
            range=min(range, N);
            for (int i=N-range; i<=N+range; i++) {
                outfile<<setprecision(6)<<setw(10)<<chop(this->Eigen_Sets[i].eval)<<",   "<<printparity(this->Eigen_Sets[i].parity);
                if (printtrans) {
                    outfile<<",   "<<setw(6)<<chop(this->Eigen_Sets[i].trans);
                    outfile<<",   "<<setprecision(5)<<setw(10)<<chop(this->Eigen_Sets[i].singletrans1);
                    outfile<<",   "<<setprecision(5)<<setw(10)<<chop(this->Eigen_Sets[i].singletrans2);
                }
                if (printes) {
                    outfile<<",   "<<setw(15)<<chop(this->Eigen_Sets[i].entanglement);
                }
                outfile<<endl;
            }
        }
        else {
            N=N/2;
            range=min(range-1, N-1);
            for (int i=N-range-1; i<=N+range; i++) {
                outfile<<setprecision(6)<<setw(15)<<chop(this->Eigen_Sets[i].eval)<<",   "<<printparity(this->Eigen_Sets[i].parity);
                if (printtrans) {
                    outfile<<",   "<<setw(6)<<chop(this->Eigen_Sets[i].trans);
                    outfile<<",   "<<setprecision(5)<<setw(10)<<chop(this->Eigen_Sets[i].singletrans1);
                    outfile<<",   "<<setprecision(5)<<setw(10)<<chop(this->Eigen_Sets[i].singletrans2);
                }
                if (printes) {
                    outfile<<",   "<<setw(15)<<chop(this->Eigen_Sets[i].entanglement);
                }
                outfile<<endl;
            }
        }
    }
    else {
        for (int i=0; i<this->Eigen_Sets.size(); i++) {
            outfile<<setprecision(6)<<setw(10)<<chop(this->Eigen_Sets[i].eval)<<",   "<<printparity(this->Eigen_Sets[i].parity);
            if (printtrans) {
                outfile<<",   "<<setw(6)<<chop(this->Eigen_Sets[i].trans);
                outfile<<",   "<<setprecision(5)<<setw(10)<<chop(this->Eigen_Sets[i].singletrans1);
                outfile<<",   "<<setprecision(5)<<setw(10)<<chop(this->Eigen_Sets[i].singletrans2);
            }
            if (printes) {
                outfile<<",   "<<setw(15)<<chop(this->Eigen_Sets[i].entanglement);
            }
            outfile<<endl;
        }
    }
    
    outfile.close();
}
void Scars::Print_Eigen_Sets(int range, string filename, bool printes, bool printph, bool printtrans){
    ofstream outfile(filename+".txt");
    ofstream outfile1p(filename+"_trans1p.txt");
    ofstream outfile1m(filename+"_trans1m.txt");
    ofstream outfile2(filename+"_trans2p.txt");
    
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
            if (abs(1.-this->Eigen_Sets[i].singletrans1)<1e-5) {
                outfile1p<<setprecision(6)<<setw(15)<<chop(this->Eigen_Sets[i].eval)<<"   "<<this->Eigen_Sets[i].parity;
                if (printph) {
                    outfile1p<<"   "<<setw(6)<<chop(this->Eigen_Sets[i].particlehole);
                }
                if (printes) {
                    outfile1p<<"   "<<setw(15)<<chop(this->Eigen_Sets[i].entanglement);
                }
                outfile1p<<endl;
            }
            if (abs(1.+this->Eigen_Sets[i].singletrans1)<1e-5) {
                outfile1m<<setprecision(6)<<setw(15)<<chop(this->Eigen_Sets[i].eval)<<"   "<<this->Eigen_Sets[i].parity;
                if (printph) {
                    outfile1m<<"   "<<setw(6)<<chop(this->Eigen_Sets[i].particlehole);
                }
                if (printes) {
                    outfile1m<<"   "<<setw(15)<<chop(this->Eigen_Sets[i].entanglement);
                }
                outfile1m<<endl;
            }
            if (abs(1.-this->Eigen_Sets[i].singletrans2)<1e-5) {
                outfile2<<setprecision(6)<<setw(15)<<chop(this->Eigen_Sets[i].eval)<<"   "<<this->Eigen_Sets[i].parity;
                if (printph) {
                    outfile2<<"   "<<setw(6)<<chop(this->Eigen_Sets[i].particlehole);
                }
                if (printes) {
                    outfile2<<"   "<<setw(15)<<chop(this->Eigen_Sets[i].entanglement);
                }
                outfile2<<endl;
            }
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
            if (abs(1.-this->Eigen_Sets[i].singletrans1)<1e-5) {
                outfile1p<<setprecision(6)<<setw(15)<<chop(this->Eigen_Sets[i].eval)<<"   "<<this->Eigen_Sets[i].parity;
                if (printph) {
                    outfile1p<<"   "<<setw(6)<<chop(this->Eigen_Sets[i].particlehole);
                }
                if (printes) {
                    outfile1p<<"   "<<setw(15)<<chop(this->Eigen_Sets[i].entanglement);
                }
                outfile1p<<endl;
            }
            if (abs(1.+this->Eigen_Sets[i].singletrans1)<1e-5) {
                outfile1m<<setprecision(6)<<setw(15)<<chop(this->Eigen_Sets[i].eval)<<"   "<<this->Eigen_Sets[i].parity;
                if (printph) {
                    outfile1m<<"   "<<setw(6)<<chop(this->Eigen_Sets[i].particlehole);
                }
                if (printes) {
                    outfile1m<<"   "<<setw(15)<<chop(this->Eigen_Sets[i].entanglement);
                }
                outfile1m<<endl;
            }
            if (abs(1.-this->Eigen_Sets[i].singletrans2)<1e-5) {
                outfile2<<setprecision(6)<<setw(15)<<chop(this->Eigen_Sets[i].eval)<<"   "<<this->Eigen_Sets[i].parity;
                if (printph) {
                    outfile2<<"   "<<setw(6)<<chop(this->Eigen_Sets[i].particlehole);
                }
                if (printes) {
                    outfile2<<"   "<<setw(15)<<chop(this->Eigen_Sets[i].entanglement);
                }
                outfile2<<endl;
            }
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
        
        //cycle_bit(statein, stateout, this->No, 1, sign);
        stateout=cycle_bits(statein, this->No);        

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
    //cout<<"translation matrix=\n"<<this->Transmat<<endl;
    
    if (this->No%2==0) {
        this->Transmat_halfL=Trans_pow(this->No/2);
    }
}
bool Scars::make_invtrans_proj(bool inv, bool trans, bool zerosector){
    Eigen::MatrixXd InvTransMat, invmatrix;
    
//    cout<<"transmat.cols,row="<<this->Transmat.cols()<<" "<<this->Transmat.rows()<<endl;
//    cout<<"transmat_half.cols,row="<<this->Transmat_halfL.cols()<<" "<<this->Transmat_halfL.rows()<<endl;
//    exit(0);
    
    if (inv) {
        invmatrix=this->P_pls;
        if (zerosector) {
            InvTransMat=this->P_pls.adjoint()*this->Transmat*this->P_pls;
        }
        else {
            InvTransMat=this->P_pls.adjoint()*this->Transmat_halfL*this->P_pls;
        }
    }
    else if (!inv) {
        invmatrix=this->P_min;
        if (zerosector) {
            InvTransMat=this->P_min.adjoint()*this->Transmat*this->P_min;
        }
        else {
            InvTransMat=this->P_min.adjoint()*this->Transmat_halfL*this->P_min;
        }
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
        //if (!(inv ^ this->one) and InvTransEval[i]>0) {
        //    Eval_N++;
        //}
        //else if ((inv ^ this->one) and InvTransEval[i]>0) {
        //    Eval_N++;
        //}
        if ((inv and trans and InvTransEval[i]>0) or (!inv and !trans and InvTransEval[i]>0)) {
            Eval_N++;
        }
        else if ((inv and !trans and InvTransEval[i]<0) or (!inv and trans and InvTransEval[i]<0)) {
            Eval_N++;
        }
    }
    
    Eigen::MatrixXd invtrans_proj0=Eigen::MatrixXd::Zero(this->bitlist.size(), Eval_N);
    this->invtrans_proj=Eigen::SparseMatrix<double>(this->bitlist.size(), Eval_N);
    //vector<Eigen::Triplet<double>> triplets;
    int counter=0;
    for (unsigned i=0; i<InvTransEval.size(); i++) {
        //if ((inv and trans and InvTransEval[i]>0) or (!inv and !trans and InvTransEval[i]>0)) {
        //    Eval_N++;
        //}
        //else if ((inv and !trans and InvTransEval[i]<0) or (!inv and trans and InvTransEval[i]<0)) {
        //    Eval_N++;
        //}
    
        if ((inv and trans and InvTransEval[i]>0) or (inv and !trans and InvTransEval[i]<0)) {
            //for (int j=0; j<this->P_pls.rows(); j++) {
            //    triplets.push_back(Eigen::Triplet<double>)(counter++, j, (this->P_pls*InvTransEvec[i])(j));
            //}
            invtrans_proj0.col(counter++)=this->P_pls*InvTransEvec[i];
        }
        else if ((!inv and trans and InvTransEval[i]<0) and (!inv and !trans and InvTransEval[i]>0)) {
            invtrans_proj0.col(counter++)=this->P_min*InvTransEvec[i];
            //for (int j=0; j<this->P_min.rows(); j++) {
            //    triplets.push_back(Eigen::Triplet<double>)(counter++, j, (this->P_min*InvTransEvec[i])(j));
            //}
        }
    }
    //invtrans_proj.ajd() * invtrans_proj = Identity.
    
    cout<<"projdim="<<invtrans_proj.cols()<<" "<<invtrans_proj.rows()<<endl; //exit(0);
    
    this->invtrans_proj=invtrans_proj0.sparseView();
    
    if (invtrans_proj.cols()==0) {
        return false;
    }
    else {
        return true;
    }
}
void Scars::makeScarShrinker_trans(int Ky){
    //Ky=1, -1 corresponds to 0 and pi sector.
    if (Ky!=1 and Ky!=-1) {
        cout<<"makeScarShrinker, Ky=1 or -1."<<endl;
    }
    
    vector<Eigen::Triplet<double>> triplets, temptrips;
    //'col' will be the ind in the two-K list. 'index' is ind in one-K list.
    int index, col=0, phase, temp, ind;
    vector<int> temp_b(this->No, 0);
    //a list that stores the index of states (and it's permutations) found.
    vector<int> found_states; state out_b;
    vector<state_int>::iterator it;
    this->statelist_2K.clear(); this->bitlist_2K.clear();
    
    for(int i=0; i<this->statelist_M.size(); i++){
        if (this->bitlist[i]==0 and Ky==-1 and statistics=="fermion") continue;//New.
        if(find(found_states.begin(),found_states.end(),i)!=found_states.end()) continue;
        
        phase=0; temptrips.clear();
        if (statistics=="fermion") {
            it=this->bitlist.begin()+i;
            temp=this->bitlist[i];
            while(true){
                index=it-this->bitlist.begin();
                double value; if (Ky==1) value=1.; else if (phase%2==0) value=1.; else value=-1.;
                temptrips.push_back(Eigen::Triplet<double>(col,index,value));
                found_states.push_back(index);
                
                temp=cycle_bits(temp, this->No);
                
                if(temp==this->bitlist[i]) break;
                
                //it=lower_bound(this->bitlist.begin(),this->bitlist.end(),temp);
                //TODO: can "find" be improved by "lower_bound"???
                it=find(this->bitlist.begin(),this->bitlist.end(),temp);
                phase++;
                if (phase>this->No) {
                    cout<<"BUG, QUIT"<<endl;
                    exit(0);
                }
            }
            
        }
        else if (statistics=="boson") {
            cout<<"cannot do boson for the time being."<<endl;
            exit(0);
        }
        
        //some cases only exist in certain momentum sectors (e.g. 0101)
        //if( (this->Np%2!=0 && Ky%(this->Np/temptrips.size())) || (this->Np%2==0 && (Ky-this->Np/2)%(this->Np/temptrips.size()))) continue;
        //TODO: note here I switched row and col, as different from Scott's code. And I did NOT take adjoint in the end.
        //cout<<"temptrips.size()="<<temptrips.size()<<endl;
        if (temptrips.size()%2==1 and Ky==-1) {
            continue;
        }
        else {
            //insert into statelist_2K;
            statelist_2K.push_back(this->statelist_M[i]);
            bitlist_2K.push_back(this->bitlist[i]);
            
            for(unsigned int j=0;j<temptrips.size();j++) {
                temptrips[j]=Eigen::Triplet<double>(temptrips[j].row(), temptrips[j].col(), temptrips[j].value()/sqrt(temptrips.size()));
            }
            triplets.insert(triplets.end(),temptrips.begin(),temptrips.end());
            col++;
        }
        
    }
    
    this->shrinkMatrix_trans=Eigen::SparseMatrix<double>(col, this->statelist_M.size());
    //cout<<statelist_2K.size()<<" "<<col<<" "<<statelist_M.size()<<endl; //exit(0);
    this->shrinkMatrix_trans.setFromTriplets(triplets.begin(),triplets.end());
    
    //cout<<this->shrinkMatrix_trans*this->shrinkMatrix_trans.adjoint()<<endl;//this is identity. both 1K and 2K state norm to 1. shrinkmatrix is not unitary (row neq col).
    //cout<<this->shrinkMatrix_trans.adjoint()*this->shrinkMatrix_trans<<endl;//this is not identity.
    //@@@@@@In my convention, the shrinkMatrix is a (2K index * 1K index) matrix. And (shrinkMatrix)*1K-state -> 2K-state.
}
void Scars::makeScarShrinker_trans_inv(bool inv){
    this->statelist_M_trans_int=this->statelist_2K;
    this->bitlist_trans_int=this->bitlist_2K;
    this->shrinkMatrix_trans_inv=this->shrinkMatrix_trans;
    
    //return ;
    
    vector<int> avoidpoints; avoidpoints.clear();
    vector<int> avoidpoints2; avoidpoints2.clear();
    
    for (int i=0; i<this->bitlist_2K.size(); i++) {
        state_int state=invert_bits(this->bitlist_2K[i], this->No);
        vector<state_int>::iterator binaryit=this->bitlist.begin();
        binaryit=lower_bound(this->bitlist.begin(), this->bitlist.end(), state);
        //state = inverted-i; j = its index.
        int j=binaryit-this->bitlist.begin();
        
        //this->statelist_M_trans_int=this->statelist_2K;
        //this->bitlist_trans_int=this->bitlist_2K;
        
        if (this->shrinkMatrix_trans.coeff(i, j)>0 and !inv) {
            //torestore.push_back(i);
            //this->statelist_M_trans_int.push_back(this->statelist_2K[i]);
            //this->bitlist_trans_int.push_back(this->bitlist_2K[i]);
            //this->shrinkMatrix_trans_inv=SparseRemoveRow(this->shrinkMatrix_trans_inv, i, this->bitlist_trans_int, this->statelist_M_trans_int);
            avoidpoints.push_back(i);
        }
        else if (this->shrinkMatrix_trans.coeff(i, j)<0 and inv) {
            //torestore.push_back(i);
            //this->statelist_M_trans_int.push_back(this->statelist_2K[i]);
            //this->bitlist_trans_int.push_back(this->bitlist_2K[i]);
            //this->shrinkMatrix_trans_inv=SparseRemoveRow(this->shrinkMatrix_trans_inv, i, this->bitlist_trans_int, this->statelist_M_trans_int);
            avoidpoints.push_back(i);
        }
        else if (this->shrinkMatrix_trans.coeff(i, j)==0) {
            //continue ;
            //cout<<"i, j = "<<i<<" "<<j<<endl;
            //print_bit(bitlist_2K[i], No);
            //print_bit(bitlist[j], No);
            
            vector<state_int>::iterator it;
            //int ind;
            state_int temp=state;
            
            while (true) {
                //temp = translated-i.
                temp=cycle_bits(temp, this->No);
                if (temp==state) {
                    cout<<"not found"<<endl; print_bit(state, No);
                    exit(0);
                }
                else {
                    it=find(this->bitlist_2K.begin(),this->bitlist_2K.end(),temp);
                    if(it!=this->bitlist_2K.end()) {
                        int ind=it-this->bitlist_2K.begin();
                        if (find(avoidpoints2.begin(),avoidpoints2.end(),i)==avoidpoints2.end()) {
                            if (inv) {
                                this->shrinkMatrix_trans_inv=SparseLinearComb(this->shrinkMatrix_trans_inv, i, ind, +1., +1.);
                            }
                            else {
                                this->shrinkMatrix_trans_inv=SparseLinearComb(this->shrinkMatrix_trans_inv, i, ind, +1., -1.);
                            }
                            avoidpoints2.push_back(ind);
                        }
                        
                        break;
                        
                    }
                    else { }
                }
            }
        }
    }
    
    Eigen::SparseMatrix<double> shrinktemp=Eigen::SparseMatrix<double>(this->shrinkMatrix_trans.cols(), this->shrinkMatrix_trans.rows()-avoidpoints.size()-avoidpoints2.size());
    bitlist_trans_int.clear();
    statelist_M_trans_int.clear();
    
    int counter=0;
    for (int i=0; i<shrinkMatrix_trans.rows(); i++) {
        if (find(avoidpoints.begin(),avoidpoints.end(),i)==avoidpoints.end() and find(avoidpoints2.begin(),avoidpoints2.end(),i)==avoidpoints2.end()) {
            shrinktemp.col(counter++)=this->shrinkMatrix_trans_inv.transpose().col(i);
            bitlist_trans_int.push_back(bitlist_2K[i]);
            statelist_M_trans_int.push_back(statelist_2K[i]);
        }
    }
    this->shrinkMatrix_trans_inv=shrinktemp.transpose();
    //cout<<this->shrinkMatrix_trans.cols()<<" "<<this->shrinkMatrix_trans.rows()<<endl;
    //cout<<this->shrinkMatrix_trans_inv.cols()<<" "<<this->shrinkMatrix_trans_inv.rows()<<endl;
    //cout<<"shrinmatrix =\n"<<this->shrinkMatrix_trans<<endl;
    //cout<<"shrinmatrix2=\n"<<this->shrinkMatrix_trans_inv<<endl;
    //cout<<"inv="<<inv<<"\n"<<this->shrinkMatrix_trans_inv*this->shrinkMatrix_trans_inv.adjoint()<<endl;
    //exit(0);
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
//Time Evolution.
void Scars::Time_Evolution(const double dt, const int Nt, const state_int& state){
    cout<<"Time Evolution is for Z2 state 1010101... for the time being."<<endl;
    
    vector<state_int>::iterator binaryit=this->bitlist.begin();
    Eigen::VectorXcd weight=Eigen::VectorXcd::Zero(this->All_States.size());
    binaryit=lower_bound(this->bitlist.begin(), this->bitlist.end(), state);
    if (binaryit!=this->bitlist.end() and *binaryit==state) {
        //cout<<"ind = "<<binaryit-this->bitlist.begin()<<endl;
        for (int i=0; i<this->All_States.size(); i++) {
            weight(i)=this->All_States[i](binaryit-this->bitlist.begin());
        }
    }
    else {
        cout<<"cannot find input state in time evolution"<<endl;
        exit(0);
    }
    //cout<<"weight = \n"<<weight<<endl;
    
    //entanglement.
    vector<vector<double>> entanglement;
    this->ee_setup(0, this->No/2, this->bitlist);
    
    int nthreads=omp_get_max_threads();
    vector<vector<vector<double>>> entanglement_t(nthreads);
    
    #pragma omp parallel for
    for (int i=0; i<Nt; i++) {
        int coren=omp_get_thread_num();
        double t=dt*i;
        
        //calculate ee.
        Eigen::VectorXcd state_t=Eigen::VectorXcd::Zero(this->bitlist.size());
        for (int k=0; k<this->bitlist.size(); k++) {
            state_t+=weight(k)*polar(1., -1.*this->Eval[k]*t)*this->All_States[k];
        }
        Eigen::MatrixXcd rho2;
        this->ee_compute_rho(state_t, rho2, this->bitlist);
        //entanglement=this->ee_eval_rho(rho2);
        
        //cout<<"\nweight_t = \n"<<Eigen_Mcd_chop(weight_t)<<endl;
        //cout<<"t = "<<setw(3)<<t<<" entanglement = "<<entanglement<<endl;
        //cout<<"norm = "<<weight_t.norm()<<endl;
        
        entanglement_t[coren].push_back(vector<double>{t, this->ee_eval_rho(rho2)});
    }
    for (int i=0; i<nthreads; i++) {
        for (int j=0; j<entanglement_t[i].size(); j++) {
            entanglement.push_back(entanglement_t[i][j]);
        }
    }
    sort(entanglement.begin(), entanglement.end(), compare_vecdouble);
    
    for (int i=0; i<entanglement.size(); i++) {
        cout<<"t = "<<setw(3)<<entanglement[i][0]<<" entanglement = "<<entanglement[i][1]<<endl;
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
    Eigen::MatrixXd MPS_Pmat_01=Eigen::MatrixXd::Zero(2,2);
    MPS_Pmat_01<<P00, P01, P10, P11;
    
    
    
    cout<<"mps.run()"<<endl;
    cout<<"mps0*mps1="<<chop(mpswf0.dot(mpswf1))<<endl;
    
    cout<<"\n\n@@@ Info of MPS0"<<endl;
    cout<<"<P>="<<chop(this->mpswf0.dot(this->Invmat*this->mpswf0))<<endl;
    cout<<"<H>="<<chop(this->mpswf0.dot(this->H_spamatrix*this->mpswf0))<<endl;
    cout<<"<H^2>-<H>^2="<<chop(this->mpswf0.dot(this->H_spamatrix*this->H_spamatrix*this->mpswf0)-pow(chop(this->mpswf0.dot(this->H_spamatrix*this->mpswf0)),2))<<endl;
    cout<<"entanglement spectrum = \n"<<this->get_espectrum(this->mpswf0, this->bitlist)<<endl;
    cout<<"entropy="<<this->get_entropy(this->mpswf0, this->bitlist)<<endl;
    
    cout<<"\n\n@@@ Info of MPS1"<<endl;
    cout<<"<P>="<<chop(this->mpswf1.dot(this->Invmat*this->mpswf1))<<endl;
    cout<<"<H>="<<chop(this->mpswf1.dot(this->H_spamatrix*this->mpswf1))<<endl;
    cout<<"<H^2>-<H>^2="<<chop(this->mpswf1.dot(this->H_spamatrix*this->H_spamatrix*this->mpswf1)-pow(chop(this->mpswf1.dot(this->H_spamatrix*this->mpswf1)),2))<<endl;
    cout<<"entanglement spectrum = \n"<<this->get_espectrum(this->mpswf1, this->bitlist)<<endl;
    cout<<"entropy="<<this->get_entropy(this->mpswf1, this->bitlist)<<endl;
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
void Scars::print_bitlist_trans_int(){
    cout<<"@@ print bitlist trans int @@"<<endl;
    for (int i=0; i<this->bitlist_trans_int.size(); i++) {
        cout<<"i = "<<setw(3)<<i<<"     ";
        print_bit(this->bitlist_trans_int[i], this->No);
    }
    cout<<"@@ done bitlist trans int @@"<<endl;
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
void quantum_scar_PXPPBCs(int no, int replu_range, int connectX){
    bool constrain=true, quickmode=true;
    int rangeS=1;
    string bound_cond="PBC", type="PXP";
    Scars qscars(no, rangeS, type, constrain, bound_cond, replu_range, connectX, quickmode);
    
    qscars.ee_setup(0, qscars.No/2, qscars.bitlist);
    qscars.generate_svdindmat();
    
    bool calculatetrans=true, calculatees=true, calculateph=true, calculatet=false, calculatett=false, calculatethalf=false, showes=true, showph=false, showtrans=false;
    //qscars.diag_scar_H_inv("both", calculatetrans, calculatees, calculateph);
    //qscars.Show_Eigen_Sets(50000, showtrans, showes, showph);
    
    string outfilename=type+"_PBC_"+to_string((long long int)(no))+"_reprg"+to_string((long long int)(replu_range))+"_Xs"+to_string((long long int)(connectX));
    if (qscars.diag_scar(true, +1, "t", calculatees, calculatet, calculatett, calculatethalf)) {
        qscars.Print_Eigen_Sets2(-1, outfilename+"_inv1_t1.txt", showes, showtrans);
    }
    if (qscars.diag_scar(false, +1, "t", calculatees, calculatet, calculatett, calculatethalf)) {
        qscars.Print_Eigen_Sets2(-1, outfilename+"_inv0_t1.txt", showes, showtrans);
    }
    if (qscars.diag_scar(true, -1, "t", calculatees, calculatet, calculatett, calculatethalf)) {
        qscars.Print_Eigen_Sets2(-1, outfilename+"_inv1_t0.txt", showes, showtrans);
    }
    if (qscars.diag_scar(false, -1, "t", calculatees, calculatet, calculatett, calculatethalf)) {
        qscars.Print_Eigen_Sets2(-1, outfilename+"_inv0_t0.txt", showes, showtrans);
    }
}
void quantum_scar_new(string filename){
    int no, rangeS;
    int replu_range, connectX, Ky;
    bool constrain, inv;
    string type, bound_cond;
    
    ifstream infile(filename);
    infile>>no>>rangeS;
    infile>>type>>constrain>>bound_cond;
    infile>>replu_range>>connectX;
    infile>>inv>>Ky;
    infile.close();
    
    quantum_scar_PXPPBCs(no, replu_range, connectX);
}
void quantum_scar(string filename){
    int no, rangeS;
    int replu_range, connectX, Ky;
    bool constrain, inv;
    string type, bound_cond;
    
    ifstream infile(filename);
    infile>>no>>rangeS;
    infile>>type>>constrain>>bound_cond;
    infile>>replu_range>>connectX;
    infile>>inv>>Ky;
    infile.close();
    
    //cout<<"no, rangeS, type, boundary_cond="<<no<<" "<<rangeS<<" "<<type<<" "<<bound_cond<<endl;
    
    Scars qscars(no, rangeS, type, constrain, bound_cond, replu_range, connectX);
    
//    qscars.makeScarShrinker_trans(1);
//    for (int i=0; i<qscars.bitlist_2K.size(); i++) {
//        print_bit(qscars.bitlist_2K[i], qscars.No);
//    }
//    cout<<"shrinkermatrix=\n"<<qscars.shrinkMatrix_trans<<endl;
//    cout<<qscars.shrinkMatrix_trans*qscars.shrinkMatrix_trans.adjoint()<<endl;//this is identity. both 1K and 2K state norm to 1. shrinkmatrix is not unitary (row neq col).
//    qscars.makeScarShrinker_trans_inv(true);
//    qscars.makeScarShrinker_trans_inv(false);
//    
//    qscars.makeScarShrinker_trans(-1);
//    for (int i=0; i<qscars.bitlist_2K.size(); i++) {
//        print_bit(qscars.bitlist_2K[i], qscars.No);
//    }
//    cout<<"shrinkermatrix=\n"<<qscars.shrinkMatrix_trans<<endl;
//    cout<<qscars.shrinkMatrix_trans*qscars.shrinkMatrix_trans.adjoint()<<endl;//this is identity. both 1K and 2K state norm to 1. shrinkmatrix is not unitary (row neq col).
//    qscars.makeScarShrinker_trans_inv(true);
//    qscars.makeScarShrinker_trans_inv(false);
//    exit(0);
    
    qscars.diag_scar_H();
    qscars.show_scar_energy(100);

//    state_int ztwostate=0;
//    int ind=0;
//    while (ind<no) {
//        ztwostate+=1 << ind;
//        ind+=2;
//    }
//
//    ztwostate=0;
//    ztwostate+=1 << 3;
//    ztwostate+=1 << 7;
    
//    print_bit(ztwostate, no);
//    //qscars.print_statelist();
//    qscars.Time_Evolution(0.1, 100, ztwostate);
//    exit(0);
    
    cout<<"$ project H into inverse sector"<<endl;
    qscars.proj_H_invsector();
    
    cout<<"H.dim()="<<qscars.bitlist.size()<<endl;
    
    //qscars.print_commutators();
    
    qscars.ee_setup(0, qscars.No/2, qscars.bitlist);
    
    bool calculatetrans=true, calculatees=true, calculateph=true, showes=true, showph=false, showtrans=true;
    qscars.diag_scar_H_inv("both", calculatetrans, calculatees, calculateph);
    qscars.Show_Eigen_Sets(50000, showtrans, showes, showph);
    
    cout<<"\nFull Symmetry, true, true"<<endl;
    bool usetrans;
    if (bound_cond=="PBC") {
        usetrans=true;
    }
    else {
        usetrans=false;
    }
    
    
    
    if (qscars.diag_scar_H_inv_halft(inv, usetrans, calculatees)) {
        qscars.Show_Eigen_Sets(50000, showtrans, showes, showph);
        //qscars.Print_Eigen_Sets(50000, bound_cond+"/E_S_Plot_16", showes, showph, true);
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
    
    
    
//    double E1=1.0, E2=2.0, S1=0.6, S2=1.0;
//    cout<<"\nPick States, E1,E2,S1,S2="<<E1<<" "<<E2<<" "<<S1<<" "<<S2<<endl;
//    vector<eigen_set> eigensets; vector<Eigen::VectorXd> spectrum;
//    qscars.pick_Eigen_Sets(E1, E2, S1, S2, eigensets, spectrum);
//    if (spectrum.size()!=1) {
//        cout<<"spectrum.size()!=1, ="<<spectrum.size()<<endl;
//        //exit(0);
//    }
//    for (int i=0; i<eigensets.size(); i++) {
//        cout<<"i="<<i<<" eng="<<eigensets[i].eval<<" ent="<<eigensets[i].entanglement<<endl;
//        cout<<"ED spectrum=\n"<<qscars.get_espectrum(eigensets[i].evec, qscars.bitlist)<<endl;
//    }
//    cout<<endl;
    
    
    
//    for (int i=0; i<spectrum.size(); i++) {
//
//        cout<<"i="<<i<<endl;
//        cout<<"ED spectrum=\n"<<this->get_espectrum(this->mpswf1, this->bitlist);
//
//
//
//        ofstream outfilees("OBC/scar_spectrum_"+to_string((long long int)(no))+"_"+to_string((long long int)(i))+".txt");
//        for (int j=0; j<spectrum[i].size(); j++) {
//            outfilees<<spectrum[i](j)<<endl;
//        }
//        outfilees.close();
//    }
    
    //qscars.mps_run();
    
    
    
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
    int no, rangeS, Ky;
    bool constrain, inv;
    string type, bound_cond;
    
    ifstream infile(filename);
    infile>>no>>rangeS;
    infile>>type>>constrain>>bound_cond;
    infile>>inv>>Ky;
    infile.close();
    
    Scars qscars(no, rangeS, type, constrain, bound_cond);
    
    qscars.ee_setup(0, qscars.No/2, qscars.bitlist);
    
    qscars.mps_run();
}
void plot_H_dim(int range, int N, string bouncond){
    ofstream outfile("dim_"+bouncond+to_string((long long int)range)+".txt");
    for (int i=4; i<N; i+=1) {
        Scars qscarp(i, 1, "PXP", true, bouncond, range, 1);
        outfile<<setw(5)<<i<<" "<<setw(5)<<qscarp.get_H_dim()<<endl;
    }
    
    outfile.close();
}
