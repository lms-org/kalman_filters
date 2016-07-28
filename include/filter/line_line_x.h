#ifndef LINE_X_H
#define LINE_X_H
#include "sgd.h"
#include "adam.h"
#include <eigen3/Eigen/Eigen>
#include <algorithm>
#include <utils.h>
#include <limits>
#include <iostream>
#include <cmath>


//typedef LineX<sgd::Adam> LineAdam;
//typedef LineX<sgd::SGDVanilla> LineSGDVan;

struct LineX:public sgd::Adam{

    bool fixX = false;
    bool fixY = false;

    float lineLength = 1;
    /**
     * @brief state first value
     * x0,y0,phi0,...phin values
     */
    void init(int numerOfSegments){
        state.resize(numerOfSegments+2);
        state.setZero();
        sg = [this](const Eigen::VectorXd &data,Eigen::MatrixXd& alpha){
            //int closestPoint; //in xy
            //double mindistance = minimum_distance(toXY(state,lineLength),data,closestPoint);
            alpha *= 1e-3;
            alpha(0,0) = 1e-2;
            alpha(1,1) = 1e-2;
            //das bringt sehr wenig
            Eigen::VectorXd derv = deriveDistance(data);
            for(int i = 2; i< derv.rows(); i++){
                //alpha(i,i) =alpha(i,i)*(1.0/std::pow((std::abs(i-closestPoint+2)+1),0.5)); //+2 um von xy in state Darstellung zu kommen
                //alpha(i,i) =alpha(i,i) *1.0/(std::pow(derv.rows()-i,1)/lineLength);// std::pow(i,0.5);
            }
            return derv;
        };
    }

    virtual void update(const Eigen::MatrixXd& data) override{
        if(data.rows()*data.cols() == 0){
            return;
        }
        double oldX = 0,oldY = 0;
        if(fixX)
            oldX = state(0);
        if(fixY)
            oldY = state(1);
        this->sgd::Adam::update(data);
        if(fixX)
            state(0) = oldX;
        if(fixY)
            state(1) = oldY;
    }

    Eigen::VectorXd deriveDistance(const Eigen::Vector2d& p){
        std::function<double(Eigen::VectorXd)> dist =[this,p](const Eigen::VectorXd x){
            //get the min distance
            Eigen::Matrix<double,Eigen::Dynamic,2> xy = toXY(x,lineLength);
            int unused;
            double minDistance = LineX::minimum_distance(xy,p,unused);
            return minDistance;
        };

        Eigen::VectorXd deriv = KALMAN_FILTER::differentiate<Eigen::VectorXd,double>(dist,state,KALMAN_FILTER::DOUBLE_H);
        return deriv;
    }


    void translate(float dx, float dy, float phi){
        Eigen::Matrix<double,Eigen::Dynamic,2> sxy = toXY(state,lineLength);
        Eigen::Matrix<double,3,3> rt;
        rt.setZero();
        rt(0,0) = std::cos(phi);
        rt(1,0) = std::sin(phi);
        rt(0,1) = -std::sin(phi);
        rt(1,1) = std::cos(phi);
        rt(0,2) = dx;
        rt(1,2) = dy;
        rt(1,2) = 1;
        //TODO
    }

    Eigen::Matrix<double,Eigen::Dynamic,2> toXY(){
        return toXY(state,lineLength);
    }


    static Eigen::Matrix<double,Eigen::Dynamic,2>  toXY(Eigen::VectorXd statePhi,const double lineLength){
        //std::cout<<"lineLength: "<<lineLength<<std::endl;
        Eigen::Matrix<double,Eigen::Dynamic,2>  stateXY(((int)statePhi.cols()*statePhi.rows())-1,2);
        stateXY(0,0) = statePhi(0);
        stateXY(0,1) = statePhi(1);
        double angle = 0;
        Eigen::Rotation2D<double> r(angle);
        for(int i = 2; i < statePhi.rows()*statePhi.cols();i++){
            r.angle() = angle;
            Eigen::Vector2d dXY(lineLength*cos(statePhi(i)),lineLength*sin(statePhi(i)));
            dXY = r*dXY;
            stateXY(i-1,0) = stateXY(i-2,0) + dXY(0);
            stateXY(i-1,1) = stateXY(i-2, 1) + dXY(1);
            angle += statePhi(i);
        }
        return stateXY;
    }

    static Eigen::VectorXd fromXY(const Eigen::VectorXd& stateXY){
        Eigen::VectorXd statePhi;
        statePhi(0) = stateXY(0);
        statePhi(1) = stateXY(1);
        for(int i = 2; i < stateXY.rows()*stateXY.cols();i+=2){
            statePhi(i) = stateXY((i-1)*2) + cos(stateXY(i));
            statePhi(i*2+1) = stateXY((i-1)*2 +1) + sin(stateXY(i)); //TODO
        }
        return statePhi;
    }

    static double minimum_distance(const Eigen::Matrix<double,Eigen::Dynamic,2>& xy, const Eigen::Vector2d& p, int& index) {
        double minDistance = std::numeric_limits<double>::max();
        for(int i = 0; i < ((int)xy.rows())-1; i++){
            double cDistance = minimum_distance(Eigen::Vector2d(xy(i,0),xy(i,1)),Eigen::Vector2d(xy(i+1,0),xy(i+1,1)),p);
            if(cDistance < minDistance){
                minDistance = cDistance;
                index = i;
            }

        }
        return minDistance;
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
        //std::cout<<"VALUES"<<std::endl;
        //std::cout<<v<<std::endl<<w<<std::endl<<p<<std::endl;
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

};


#endif //LINE_X_H
