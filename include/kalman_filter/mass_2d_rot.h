#pragma once
#include <Eigen/Eigen>
#include <kalman/ExtendedKalmanFilter.hpp>
#include <kalman/UnscentedKalmanFilter.hpp>
#include <iostream>

namespace KALMAN_FILTERS{
namespace MASS_2D_ROT_EKF{
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
M h(const S& x) const
{
    M measurement;
    //TODO what is measurement.x()???
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
        //create State Matrix
        //std::cout<<"BEGIN res: "<<res.x()<<" "<<res.y()<<std::endl;
        Eigen::Matrix3f stateMat;
        stateMat(0,0) = std::cos(x.phi());
        stateMat(0,1) = -std::sin(x.phi());
        stateMat(1,0) = std::sin(x.phi());
        stateMat(1,1) = std::cos(x.phi());
        stateMat(0,2) = x.x();
        stateMat(1,2) = x.y();
        stateMat(2,0) = 0;
        stateMat(2,1) = 0;
        stateMat(2,2) = 1;
        Eigen::Matrix3f transRot;
        //TODO dx,dy change if dAngle is not zero
        float dAngle = x.omega()*u.dt();
        float dx = x.vx()*u.dt();
        float dy = x.vy()*u.dt();
        transRot(0,0) = std::cos(dAngle);
        transRot(0,1) = -std::sin(dAngle);
        transRot(1,0) = std::sin(dAngle);
        transRot(1,1) = std::cos(dAngle);
        transRot(0,2) = dx;
        transRot(1,2) = dy;
        transRot(2,0) = 0;
        transRot(2,1) = 0;
        transRot(2,2) = 1;
        Eigen::Vector3f pos{0,0,1};
        stateMat = stateMat*transRot;
        pos = stateMat*pos;
        res.x() = pos(0);
        res.y() = pos(1);
        res.phi() = atan2(stateMat(1,0),stateMat(0,0)); //TODO this seems to be wrong but I don't know why
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

    void setMeasurementVec(const T vx, const T vy,const T omega){
        z.setZero();
        z.vx() = vx;
        z.vy() = vy;
        z.omega() = omega;
        //z.x() = 1; //TODO why does this compile, there is no method called x() is MyMeasurement?
        //z.y() = 2; //It seems to be a State???
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
