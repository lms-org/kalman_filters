#ifndef KALMAN_FILTERS_CTRA_EKF_H
#define KALMAN_FILTERS_CTRA_EKF_H
#include <eigen3/Eigen/Eigen>
#include <kalman/ExtendedKalmanFilter.hpp>
#include <kalman_filter/ctra_measurement_model.h>
#include <kalman_filter/ctra_system_model.h>

//TODO not done yet

namespace KALMAN_FILTERS{
namespace CTRA{

struct CTRAContainer{
    typedef float T;

    typedef CTRA::State<T> State;
    typedef CTRA::Control<T> Control;
    typedef CTRA::Measurement<T> Measurement;

    typedef CTRA::MeasurementModel<T> MeasurementModel;
    typedef CTRA::SystemModel<T> SystemModel;
    typedef Kalman::ExtendedKalmanFilter<State> Filter;

    Control u;
    Measurement z;

    SystemModel sys;
    MeasurementModel mm;

    Filter filter;

    void setMeasurementVec(const T v, const T ax, const T ay, const T omega){
        z.v() = v;
        z.ax() = ax;
        z.ay() = ay;
        z.omega() = omega;
    }

    void setMeasurementCov(const T axVar,const T ayVar,const T vVar,const T omegaVar){
        Kalman::Covariance< Measurement > cov;
        cov.setZero();
        cov(Measurement::AX,    Measurement::AX)    = axVar;
        cov(Measurement::AY,    Measurement::AY)    = ayVar;
        cov(Measurement::V,     Measurement::V)     = vVar;
        cov(Measurement::OMEGA, Measurement::OMEGA) = omegaVar;
        mm.setCovariance(cov);
    }

    void predict(T dt){
        u.dt = dt;
        filter.predict(sys, u);
    }

    void update(){
        filter.update(mm, z);
    }

    void pu(T dt){
        predict(dt);
        update();
    }

    //TODO was ist mit initial state covariance
    //TODO was it mit process noise covariance
};

} //namespace CTRA
} //namespace KALMAN_FILTERS_CTRA_EKF
#endif //KALMAN_FILTERS_CTRA_EKF
