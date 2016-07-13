#ifndef KALMAN_FILTERS_LINE_EKF
#define KALMAN_FILTERS_LINE_EKF
#include <eigen3/Eigen/Eigen>
#include <kalman/ExtendedKalmanFilter.hpp>

//TODO not done yet

namespace KALMAN_FILTERS{
namespace lineKalman{
typedef float T;
//State
template<typename T>
class State : public Kalman::Vector<T, Eigen::Dynamic>{
    //erster x-Wert immer gleich 0
    //-> erster wert y variabel
    //-> winkel
    static constexpr size_t Y0 = 0;
public:
    KALMAN_VECTOR(State, T, Eigen::Dynamic)
    float partLength;
    T x0(){
        return 0;
    }
    T y0(){
        return (*this)[Y0];
    }
    Kalman::Vector<T, Eigen::Dynamic> xy(){
        Kalman::Vector<T, Eigen::Dynamic> res;
        //TODO
        return res;
    }

    Kalman::Vector<T, Eigen::Dynamic> angles(){
        Kalman::Vector<T, Eigen::Dynamic> res;
        //TODO set values (1-end)
        return res;
    }
};

template<typename T>
class Control : public Kalman::Vector<T, Eigen::Dynamic>{
public:
    //hier übergibt man dx, dy, dphi?
    KALMAN_VECTOR(Control, T, Eigen::Dynamic)
};

template<typename T>
class Measurement : public Kalman::Vector<T, Eigen::Dynamic>{
    //hier üebrgibt man die gefundenen neuen Punkte
public:
    KALMAN_VECTOR(Measurement, T, Eigen::Dynamic)
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
    measurement = x;
    return measurement;
}
};
/**
 * @brief System model CTRA motion
 *
 * This is the system model defining how our car moves from one
 * time-step to the next, i.e. how the system state evolves over time.
 *
 * @param T Numeric scalar type
 * @param CovarianceBase Class template to determine the covariance representation
 *                       (as covariance matrix (StandardBase) or as lower-triangular
 *                       coveriace square root (SquareRootBase))
 */
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
        //TODO
        return x;
    }
};
class ExtendedKalmanFilterForLine{

public:
    typedef float T;

    typedef lineKalman::State<T> State;
    typedef lineKalman::Control<T> Control;
    typedef lineKalman::Measurement<T> Measurement;

    typedef lineKalman::MeasurementModel<T> MeasurementModel;
    typedef lineKalman::SystemModel<T> SystemModel;
    typedef Kalman::ExtendedKalmanFilter<State> Filter;

protected:
    float partLength;
    //Control u;
    //Measurement z;

    //SystemModel sys;
    //MeasurementModel mm;

    //Filter filter;
    void computeMeasurement();
    void computeFilterStep();
    void initFilter();

public:
    void init(const Eigen::VectorXf &points, const float partLength);
    void update(const Eigen::VectorXf &points,const float dx,const float dy,const float dphi);

};


} //namespace lineKalman
} //namespace KALMAN_FILTERS
#endif //KALMAN_FILTERS_LINE_EKF
