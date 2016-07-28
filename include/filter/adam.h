#ifndef SGD_ADAM_H
#define SGD_ADAM_H
#include <Eigen/Eigen>
#include "sgd.h"
#include <iostream>

namespace sgd{
struct Adam:public SGDContainer{
    /**
     * @brief m first momentum vector
     */
    Eigen::VectorXd m;
    /**
     * @brief v second momentum vector
     */
    Eigen::VectorXd v;
    /**
     * @brief b1,b2 [0,1] Exponential decay rates for the moment estimates
     */
    double b1=0.9,b2=0.999;
    /**
     * @brief a stepsize
     */
    double a=0.001;
    /**
     * @brief e
     */
    double e= 1e-8;


    /**
     * @brief update
     * @param data are col vectors (each row is one vector)
     */
    virtual void update(const Eigen::MatrixXd& data) override{
        Adam::opt(sg,state,data,numerOfIterations,m,v,b1,b2,a,e);
    }

    static void opt(SG sg,Eigen::VectorXd& state,const Eigen::MatrixXd& data, const int numberOfIterations,Eigen::VectorXd &m,Eigen::VectorXd &v, const double b1,const double b2,double a,double e){
        for(int i = 0; i < numberOfIterations; i++){
            int j = 0;
            do{
                Eigen::MatrixXd alpha = Eigen::MatrixXd::Identity(state.rows(),state.rows());
                Eigen::VectorXd derv;
                if(data.cols() <= 0)
                    derv = sg(Eigen::VectorXd(),alpha);
                else{
                    derv = sg(data.col(j),alpha);
                }
                //calculate momentums
                m = b1*m+(1-b1)*derv;
                v = b2*v+((1-b2)*derv.array()*derv.array()).matrix();

                //calculate bias corrected momentums
                Eigen::VectorXd mC = m/(1-b1);
                Eigen::VectorXd mV = v/(1-b2*b2);
                state = state - (a*mC.array()/(mV.array().sqrt()+e)).matrix();
                j++;
            }while(j < data.cols());
        }
    }
};
}//namespace sgd

#endif //SGD_ADAM_H
