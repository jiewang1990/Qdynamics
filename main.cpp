//
//  main.cpp
//  quantumscar
//
//  Created by Jie Wang on 12/10/18.
//  Copyright © 2018 Jie Wang. All rights reserved.
//

#include "quantumscar.h"
void readfile(string filename, vector<double>& energy){
    energy.clear();
    string path="/Users/jiewang/Google\ Drive/jie_programs/quantumscar/data/";
    ifstream infile(path+filename);
    string line;
    getline(infile, line);
    for (line; getline(infile, line); ) {
        size_t pos = line.find(",");
        energy.push_back(stod(line.substr(0, pos)));
    }
    infile.close();
}
vector<double> levelsta(const vector<double> energy){
    vector<double> etmp=energy, etmp2, etmpdiff, etmpdiv;
    sort(etmp.begin(), etmp.end());
    
    //delete zero modes.
    for (int i=0; i<etmp.size(); i++) {
        if (etmp[i]!=0.) {
            etmp2.push_back(etmp[i]);
        }
    }
    
    //calculate energy difference.
    for (int i=0; i<etmp2.size()-1; i++) {
        etmpdiff.push_back(etmp2[i+1]-etmp2[i]);
    }
    
    for (int i=0; i<etmpdiff.size()-1; i++) {
        if (etmpdiff[i+1]==0 or etmpdiff[i]==0) {
            etmpdiv.push_back(0.);
        }
        else {
            etmpdiv.push_back(min(etmpdiff[i+1]/etmpdiff[i], etmpdiff[i]/etmpdiff[i+1]));
        }
    }
    
    //for (int i=0; i<etmp2.size(); i++) cout<<etmp2[i]<<endl; cout<<endl;
    //for (int i=0; i<etmpdiff.size(); i++) cout<<etmpdiff[i]<<endl; cout<<endl;
    //for (int i=0; i<etmpdiv.size(); i++) cout<<etmpdiv[i]<<endl;
    
    return etmpdiv;
}
Eigen::VectorXcd psi(double lambda, double d, double phi, double dz, double eta){
    complex<double> e=eta*sqrt(d*d-lambda*lambda/4-complex<double>(0,1)*lambda*d*sin(phi));
    Eigen::VectorXcd temp=Eigen::VectorXcd::Zero(2);
    temp(0)=polar(d, phi)+lambda/2;
    temp(1)=e;
    temp/=temp.norm();
    return temp;
}
complex<double> chern(int N, double lambda, bool flip){
    complex<double> sum=0.;
    for (int i=0; i<N; i++) {
        if (flip) {
            sum+=psi(-lambda, 1, 2*M_PI*(i+1)/N, 0, 1).dot(psi(lambda, 1, 2*M_PI*i/N, 0, 1))/(1.*N);
        }
        else {
            sum+=psi(lambda, 1, 2*M_PI*(i+1)/N, 0, 1).dot(psi(lambda, 1, 2*M_PI*i/N, 0, 1))/(1.*N);
        }
    }
    return complex<double>(0, 1)*sum;
}

int main(int argc, const char * argv[]) {
    vector<double> energy;
    readfile("PXP_PBC_18_reprg1_Xs1_inv1_t1.txt", energy);
    
    vector<double> etmp=levelsta(energy);
    ofstream outfile("tilder.txt");
    for (int i=0; i<etmp.size(); i++) {
        outfile<<etmp[i]<<endl;
    }
    outfile.close();
    
    quantum_scar_new("params");
//    quantum_scar_invt("params");
//    quantum_scar_mps("params");
    
//    int N=6;
//    plot_H_dim(1, N, "OBC");
//    plot_H_dim(2, N, "PBC");
//    plot_H_dim(3, N, "OBC");
}
