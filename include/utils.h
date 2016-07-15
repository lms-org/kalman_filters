#ifndef KALMAN_FILTER_UTILS_H
#define KALMAN_FILTER_UTILS_H
#include <functional>
#include <Eigen/Eigen>
#include <cmath>


namespace KALMAN_FILTER{
    const double DOUBLE_H = 1e-7;
    template<typename M,typename T>
    M differentiate(std::function<T(M)> f,const M& x, T h){
        M res;
        M vh;
        vh.setZero(h);
        for(int r = 0; r < x.rows(); r++){
            for(int c = 0; c < x.cols(); c++){
                vh(r,c) = h;
                res(r,c)= (f(x-vh)-f(x+vh))/(2*h);
                vh(r,c) = 0;
            }
        }
        return res;
    }

    double differentiate(std::function<double(double)> f,float x, double h = 1e-7){
        return (f(x-h)-f(x+h))/(2*h);
    }
}


#endif
