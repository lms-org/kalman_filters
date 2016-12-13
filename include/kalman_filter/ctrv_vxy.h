#pragma once
#include <Eigen/Eigen>
#include <kalman/ExtendedKalmanFilter.hpp>
#include <kalman/UnscentedKalmanFilter.hpp>
#include <iostream>

namespace kalman_filters{
namespace ctrv_vxy{
typedef float T;
//State
template<typename T>
class State : public Kalman::Vector<T, 6>{
public:
    static constexpr size_t X = 0;
    static constexpr size_t Y = 1;
    static constexpr size_t VX = 2;
    static constexpr size_t VY = 3;
    /**
     * @brief PHI rotation angle
     */
    static constexpr size_t PHI = 4;
    /**
     * @brief OMEGA angular velocity
     */
    static constexpr size_t OMEGA = 5;

    KALMAN_VECTOR(State, T, 6)
    T& x(){
        return (*this)[X];
    }
    T& y(){
        return (*this)[Y];
    }
    T& vx(){
        return (*this)[VX];
    }
    T& vy(){
        return (*this)[VY];
    }
    T& phi(){
        return (*this)[PHI];
    }
    T& omega(){
        return (*this)[OMEGA];
    }
    T x() const{
        return (*this)[X];
    }
    T y() const{
        return (*this)[Y];
    }
    T vx() const{
        return (*this)[VX];
    }
    T vy() const{
        return (*this)[VY];
    }
    T phi() const{
        return (*this)[PHI];
    }
    T omega() const{
        return (*this)[OMEGA];
    }
};

template<typename T>
class Control : public Kalman::Vector<T, 1>{
public:
    KALMAN_VECTOR(Control, T, 1)

    //! time since filter was last called
    static constexpr size_t DT = 0;

    T dt()      const { return (*this)[ DT ]; }
    T& dt()     { return (*this)[ DT ]; }
};

template<typename T>
class Measurement : public Kalman::Vector<T,3>{
public:
    static constexpr size_t VX = 0;
    static constexpr size_t VY = 1;
    static constexpr size_t OMEGA = 2;
    KALMAN_VECTOR(Measurement, T, 3)
    T vx() const{
        return (*this)[VX];
    }
    T vy() const{
        return (*this)[VY];
    }
    T omega() const{
        return (*this)[OMEGA];
    }
    T& vx(){
        return (*this)[VX];
    }
    T& vy(){
        return (*this)[VY];
    }
    T& omega(){
        return (*this)[OMEGA];
    }
};

template<typename T, template<class> class CovarianceBase = Kalman::StandardBase>
class MeasurementModel : public Kalman::LinearizedMeasurementModel<State<T>, Measurement<T>, CovarianceBase>{

public:
    //! State type shortcut definition
    typedef State<T> S;

    //! Measurement type shortcut definition
    typedef Measurement<T> M;

/**
 * @brief Definition of (possibly non-linear) measurement function
 *
 * This function maps the system state to the measurement that is expected
 * to be received from the sensor assuming the system is currently in the
 * estimated state.
 *
 * @param [in] x The system state in current time-step
 * @returns The (predicted) sensor measurement for the system state
 */
M h(const S& x) const{
    M measurement;

    //We have to rotate it from the absolute system to the current one
    /*
    Eigen::Matrix2f vRot;
    vRot(0,0) = std::cos(x.phi());
    vRot(0,1) = -std::sin(x.phi());
    vRot(1,0) = std::sin(x.phi());
    vRot(1,1) = std::cos(x.phi());
    Eigen::Vector2f vV;
    vV(0) = x.vx();
    vV(1) = x.vy();
    vV = vRot*vV;
    measurement.vx() = vV(0);
    measurement.vy() = vV(1);
    measurement.omega() = x.omega();
    std::cout<<"velocity "<<measurement.vx()<<" "<<measurement.vy()<<std::endl;
    */

    measurement.vx() = x.vx();
    measurement.vy() = x.vy();
    measurement.omega() = x.omega();

    return measurement;
}

void updateJacobians(const S& x){
    (void)x;
    this->H.setIdentity();
    this->V.setIdentity();
    //TODO has to be done, but I will give UKF a try :)
}
};

template<typename T, template<class> class CovarianceBase = Kalman::StandardBase>
class SystemModel : public Kalman::LinearizedSystemModel<State<T>, Control<T>, CovarianceBase>{
public:
    //! State type shortcut definition
    typedef State<T> S;

    //! Control type shortcut definition
    typedef Control<T> C;
    /**
     * @brief Definition of (non-linear) state transition function
     *
     * This function defines how the system state is propagated through time,
     * i.e. it defines in which state \f$\hat{x}_{k+1}\f$ is system is expected to
     * be in time-step \f$k+1\f$ given the current state \f$x_k\f$ in step \f$k\f$ and
     * the system control input \f$u\f$.
     *
     * @param [in] x The system state in current time-step
     * @param [in] u The control vector input
     * @returns The (predicted) system state in the next time-step
     */
    S f(const S& x, const C& u) const{
        S res = x;
        float dAngle = x.omega()*u.dt();
        //rotate the velocity
        Eigen::Matrix2f vRot;
        vRot(0,0) = std::cos(dAngle);
        vRot(0,1) = std::sin(dAngle);
        vRot(1,0) = -std::sin(dAngle);
        vRot(1,1) = std::cos(dAngle);
        Eigen::Vector2f vV;
        vV(0) = x.vx();
        vV(1) = x.vy();
        vV = vRot*vV;
        res.vx() = vV(0);
        res.vy() = vV(1);

        res.x() += res.vx()*u.dt();
        res.y() += res.vy()*u.dt();
        res.phi() += dAngle;
        return res;
    }


    void updateJacobians( const S& x, const C& u ){
        (void)x;
        (void)u;
        this->F.setIdentity();
        this->W.setIdentity();
        //TODO has to be done, but I will give UKF a try :)
    }
};


template<class F>
class FilterForMassModel{
    typedef float T;

    typedef State<T> MyState;
    typedef Control<T> MyControl;
    typedef Measurement<T> MyMeasurement;

    typedef MeasurementModel<T> MyMeasurementModel;
    typedef SystemModel<T> MySystemModel;
    typedef F Filter; //TODO The Filter state could differ from MyState

protected:
    MyControl u;
    MyMeasurement z;

    MySystemModel sys;
    MyMeasurementModel mm;

    Filter filter;

public:
    MyState lastState;
    /**
     * @brief setMeasurementVec
     * @param vx in local coordinates
     * @param vy in local coordinates
     * @param omega
     */
    void setMeasurementVec(const T vx, const T vy,const T omega){
        z.setZero();
        //passive rot matrix but the phi is negative -> looks like active rot matrix
        Eigen::Matrix2f vRot;
        vRot(0,0) = std::cos(lastState.phi());
        vRot(0,1) = -std::sin(lastState.phi());
        vRot(1,0) = std::sin(lastState.phi());
        vRot(1,1) = std::cos(lastState.phi());
        Eigen::Vector2f vV;
        vV(0) = vx;
        vV(1) = vy;
        vV = vRot*vV;
        z.vx() = vV(0);
        z.vy() = vV(1);
        z.omega() = omega;
    }
    void predict(T dt){
        u.dt() = dt;
        lastState = filter.predict(sys, u);
    }

    void init(){
        // Init kalman
        MyState s;
        s.setZero();
        filter.init(s);

        // Set initial state covariance
        Kalman::Covariance<MyState> stateCov;
        stateCov.setIdentity();
        //TODO set stateCov
        filter.setCovariance(stateCov);

        // Set process noise covariance (Q)
        Kalman::Covariance<MyState> ncov;
        ncov.setIdentity();
        //TODO set ncov
        sys.setCovariance(ncov);

        // Set measurement covariances (R)
        Kalman::Covariance< MyMeasurement > cov;
        cov.setIdentity();
        //TODO cov.setZero();
        mm.setCovariance(cov);
    }

    void update(){
        lastState = filter.update(mm, z);
    }

    void pu(T dt){
        predict(dt);
        update();
    }
};

typedef FilterForMassModel<Kalman::ExtendedKalmanFilter<State<T>>> MassModelEKF;
typedef FilterForMassModel<Kalman::UnscentedKalmanFilter<State<T>>> MassModelUKF;

} //namespace KALMAN_FILTERS_MASS_EKF
} //namespace KALMAN_FILTERS
