#ifndef LMS_SGD_H
#define LMS_SGD_H

#include <eigen3/Eigen/Eigen>
#include <functional>
namespace sgd{

typedef std::function<Eigen::VectorXd(const Eigen::VectorXd &data,Eigen::MatrixXd& alpha)> SG;


struct SGDContainer{
    Eigen::VectorXd state;
    int numerOfIterations = 1;
    SG sg;

    virtual void update(const Eigen::MatrixXd& data)= 0;

    virtual ~SGDContainer(){}
};

struct SGDVanilla:public SGDContainer{
    /**
     * @brief update
     * @param data are col vectors (each row is one vector)
     */
    virtual void update(const Eigen::MatrixXd& data) override{
        opt(sg,state,data,numerOfIterations);
    }

    static void opt(SG sg,Eigen::VectorXd& state,const Eigen::MatrixXd& data, const int numberOfIterations){
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
                state = state -alpha*derv;
                j++;
            }while(j < data.cols());
        }
    }
    virtual ~SGDVanilla(){}

};

}

#endif //SGD_H
