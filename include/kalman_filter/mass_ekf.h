#ifndef KALMAN_FILTERS_MASS_EKF_H
#define KALMAN_FILTERS_MASS_EKF_H
#include <eigen3/Eigen/Eigen>
#include <kalman/ExtendedKalmanFilter.hpp>

//TODO not done yet

namespace KALMAN_FILTERS{
namespace MASS_EKF{
typedef float T;
//State
template<typename T>
class State : public Kalman::Vector<T, 4>{
    //erster x-Wert immer gleich 0
    //-> erster wert y variabel
    //-> winkel
    static constexpr size_t X = 0;
    static constexpr size_t Y = 1;
    static constexpr size_t VX = 2;
    static constexpr size_t VY = 3;
public:
    KALMAN_VECTOR(State, T, 4)
    T x(){
        return (*this)[X];
    }
    T y(){
        return (*this)[Y];
    }
    T vx(){
        return (*this)[VX];
    }
    T vy(){
        return (*this)[VY];
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
class Measurement : public Kalman::Vector<T, 2>{
    //hier üebrgibt man die gefundenen neuen Punkte
public:
    static constexpr size_t X = 0;
    static constexpr size_t Y = 1;
    KALMAN_VECTOR(Measurement, T, 2)
    T x() const{
        return (*this)[X];
    }
    T y() const{
        return (*this)[Y];
    }
    T& x(){
        return (*this)[X];
    }
    T& y(){
        return (*this)[Y];
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
    //TODO wie bekomme ich dann hier dx,dy,dhi rein oder brauch ich die hier garnicht?
    measurement.x() = x.x();
    measurement.y() = x.y();
    return measurement;
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
        res.x() = x.x()+ x.vx()*u.dt();
        res.y() = x.y()+x.vy()*u.dt();
        return res;
    }


    void updateJacobians( const S& x, const C& u ){
        this->F.setIdentity();
        this->F(S::X,S::VX) = u.dt();
        this->F(S::Y,S::VY) = u.dt();
    }
};
class ExtendedKalmanFilterForMassModel{

public:
    typedef float T;

    typedef State<T> MyState;
    typedef Control<T> MyControl;
    typedef Measurement<T> MyMeasurement;

    typedef MeasurementModel<T> MyMeasurementModel;
    typedef SystemModel<T> MySystemModel;
    typedef Kalman::ExtendedKalmanFilter<MyState> Filter;

protected:
    MyControl u;
    MyMeasurement z;

    MySystemModel sys;
    MyMeasurementModel mm;

    Filter filter;

public:

    void setMeasurementVec(const T x, const T y){
        z.x() = x;
        z.y() = y;
    }
    void predict(T dt){
        u.dt() = dt;
        filter.predict(sys, u);
    }
    //TODO
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
        filter.update(mm, z);
    }

    void pu(T dt){
        predict(dt);
        update();
    }



};


} //namespace KALMAN_FILTERS_MASS_EKF
} //namespace KALMAN_FILTERS
#endif //KALMAN_FILTERS_MASS_EKF_H
