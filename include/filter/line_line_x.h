#ifndef LINE_X_H
#define LINE_X_H

#include <eigen3/Eigen/Eigen>
#include <algorithm>
#include <utils.h>
#include <limits>

struct LineX{

    float lineLength = 1;
    /**
     * @brief state first value
     * x0,y0,phi0,...phin values
     */
    Eigen::VectorXd state;

    static Eigen::VectorXd toXY(Eigen::VectorXd statePhi){
        Eigen::VectorXd stateXY;
        stateXY(0) = statePhi(0);
        stateXY(1) = statePhi(1);
        for(int i = 2; i < statePhi.rows();i++){
            stateXY(i*2) = stateXY((i-1)*2) + cos(statePhi(i));
            stateXY(i*2+1) = stateXY((i-1)*2 +1) + sin(statePhi(i));
        }
        return stateXY;
    }

    static Eigen::VectorXd fromXY(const Eigen::VectorXd& stateXY){
        Eigen::VectorXd statePhi;
        statePhi(0) = stateXY(0);
        statePhi(1) = stateXY(1);
        for(int i = 2; i < stateXY.rows()*stateXY.cols();i+=2){
            statePhi(i) = stateXY((i-1)*2) + cos(stateXY(i));
            statePhi(i*2+1) = stateXY((i-1)*2 +1) + sin(stateXY(i));
        }
        return statePhi;
    }
    /**
     * @brief minimum_distance
     * from http://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
     * @param v
     * @param w
     * @param p
     * @return
     */
    static double minimum_distance(const Eigen::Vector2d& v, const Eigen::Vector2d& w, const Eigen::Vector2d& p) {
        // Return minimum distance between line segment vw and point p
        const double l2 = (v-w).squaredNorm();//length_squared(v, w);  // i.e. |w-v|^2 -  avoid a sqrt
        if (l2 == 0.0) return (p-v).norm();   // v == w case
        // Consider the line extending the segment, parameterized as v + t (w - v).
        // We find projection of point p onto the line.
        // It falls where t = [(p-v) . (w-v)] / |w-v|^2
        // We clamp t from [0,1] to handle points outside the segment vw.
        const float t = std::max(0.0,double(std::min(1.0, (p - v).dot(w - v) / l2)));
        const Eigen::Vector2d projection = v + t * (w - v);  // Projection falls on the segment
        return (p-projection).norm();
    }

    static double distanceToLine(const Eigen::Vector2d& v, const Eigen::Vector2d& w, const Eigen::Vector2d& p){
        // Return minimum distance between line segment vw and point p
        const double l2 = (v-w).squaredNorm();//length_squared(v, w);  // i.e. |w-v|^2 -  avoid a sqrt
        if (l2 == 0.0) return (p-v).norm();   // v == w case
        // Consider the line extending the segment, parameterized as v + t (w - v).
        // We find projection of point p onto the line.
        // It falls where t = [(p-v) . (w-v)] / |w-v|^2
        const float t = (p - v).dot(w - v) / l2;
        const Eigen::Vector2d projection = v + t * (w - v);  // Projection falls on the segment
        return (p-projection).norm();
    }

    void deriveDistance(const Eigen::Vector2d& p){
        std::function<double(Eigen::VectorXd)> dist =[p](const Eigen::VectorXd x){
            //get the min distance
            double minDistance = std::numeric_limits<double>::max();
            Eigen::VectorXd xy = toXY(x);
            for(int i = 0; i < x.rows(); i +=2){
                double cd = minimum_distance(Eigen::Vector2d(xy(i),xy(i+1)),Eigen::Vector2d(xy(i+2),xy(i+3)),p);
                if(cd < minDistance){
                    minDistance = cd;
                }
            }
            return minDistance;
        };

        Eigen::VectorXd deriv = KALMAN_FILTER::differentiate<Eigen::VectorXd,double>(dist,state,KALMAN_FILTER::DOUBLE_H);
    }


    void translate(float dx, float dy, float dphi){

    }

};



#endif //LINE_X_H
