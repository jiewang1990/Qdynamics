//
//  main.cpp
//  quantumscar
//
//  Created by Jie Wang on 12/10/18.
//  Copyright © 2018 Jie Wang. All rights reserved.
//

#include "quantumscar.h"
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
    quantum_scar_new("params");
//    quantum_scar_invt("params");
//    quantum_scar_mps("params");
    
//    int N=6;
//    plot_H_dim(1, N, "OBC");
//    plot_H_dim(2, N, "PBC");
//    plot_H_dim(3, N, "OBC");
}
