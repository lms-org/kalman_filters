#include <kalman_filter/mass_ekf.h>
#include <filter/line_line_x.h>
#include <filter/sgd.h>
#include <iostream>
#include <Eigen/Eigen>
#include <cmath>


#include <QtGui/QGuiApplication>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtCharts/QChartView>
#include <QtCharts/QSplineSeries>
#include <QtCharts/QScatterSeries>
#include <filter/adam.h>

int main(int argc, char** argv){
    (void)argc;
    (void)argv;
    std::function<double(double)> f1 = [](double x){
        return x*x;
    };

    std::cout<<(KALMAN_FILTER::differentiate(f1,3))<<std::endl;
    std::cout<<"###"<<std::endl;

    std::function<double(Eigen::VectorXd)> f2 =[](Eigen::VectorXd x){
        return cos(x(0))+sin(x(1))+x(2)*x(2);
    };
    Eigen::VectorXd pos(3);
    pos(0)=M_PI/2.0;
    pos(1)=M_PI/2.0;
    pos(2)=2;
    std::cout<<(KALMAN_FILTER::differentiate<Eigen::VectorXd,double>(f2,pos,KALMAN_FILTER::DOUBLE_H))<<std::endl;


    LineX lineX;
    lineX.init(8);
    lineX.lineLength = 0.2;
    lineX.fixX = true;
    lineX.fixY = false;

    //PROBLEM: Punkte, welche länger als lineX sind ziehen es in die länge bei Knicken


    //QWT http://qwt.sourceforge.net/
    //oder einfach nur qt5 http://doc.qt.io/qt-5/qtdatavisualization-index.html

    QT_CHARTS_USE_NAMESPACE

    QApplication app(argc, argv);

    QLineSeries fit;
    QScatterSeries points;

    QChart *chart = new QChart();
    chart->legend()->hide();
    chart->addSeries(&points);
    chart->addSeries(&fit);
    chart->createDefaultAxes();
    chart->setTitle("Line X Model #Ada");

    QChartView chartView(chart);
    chartView.setRenderHint(QPainter::Antialiasing);

    QMainWindow window;
    window.setCentralWidget(&chartView);
    window.resize(400, 300);
    window.show();
    //chart->axisX()->setRange(0,10);
    //chart->axisY()->setRange(-1,-1);
    int iter = 0;
    while(true){
        iter++;
        std::cout <<"iteration "<<iter<<std::endl;
        chart->setTitle(QString::fromStdString("Line X Model Iteration: " + std::to_string(iter)));
        //create new measurement-points
        int numberOfPoints = 10;
        Eigen::Matrix<double,Eigen::Dynamic,2> mData(numberOfPoints,2);
        double d = 0.2;
        double offset = 0.2;
        float randomF = 1;

        //range
        float xmin = 0;
        float xmax = 0;
        float ymin = 0;
        float ymax = 0;
        for(int i = 0; i < numberOfPoints; i++){
            float r = ((double) rand() / (RAND_MAX));
            mData(i,0) = 0.2*i;
#define PARABEL
#ifdef KREUZUNG
            //kreuzung
            if(i < 3){
                mData(i,1) = 2;
            }else{
                mData(i,0) = 3;
                mData(i,1) = mData(i-1,1)-1;
            }
#endif
#ifdef PARABEL
            //parabel
            mData(i,1) = d*mData(i,0)*mData(i,0)+offset + randomF*r;
#endif
            //range checks
            if(mData(i,0) < xmin){
                xmin = mData(i,0);
            }
            if(mData(i,0) > xmax){
                xmax = mData(i,0);
            }
            if(mData(i,1) < ymin){
                ymin = mData(i,1);
            }
            if(mData(i,1) > ymax){
                ymax = mData(i,1);
            }
        }
        //std::cout << lineX.state<<std::cout;
        //fit line
        double error = 0;
        for(int i = 0; i < numberOfPoints; i++){
            lineX.update(Eigen::Vector2d(mData(i,0),mData(i,1)));
        }
        error = error/numberOfPoints;
        //std::cout<<"error: " <<error<<std::endl;

        //conversion to qt points
        //std::cout <<"state: "<<lineX.state<<std::endl;
        Eigen::Matrix<double,Eigen::Dynamic,2> xy = lineX.toXY();
        //std::cout <<"statexy: "<<xy<<std::endl;

        points.clear();
        fit.clear();
        for(int i = 0; i < xy.rows(); i++){
            fit.append(xy(i,0),xy(i,1));
            //range checks
            if(xy(i,0) < xmin){
                xmin = xy(i,0);
            }
            if(xy(i,0) > xmax){
                xmax = xy(i,0);
            }
            if(xy(i,1) < ymin){
                ymin = xy(i,1);
            }
            if(xy(i,1) > ymax){
                ymax = xy(i,1);
            }
        }
        for(int i = 0; i < numberOfPoints; i++){
            points.append(mData(i,0),mData(i,1));
        }

        //std::cout <<"points: "<<mData<<std::endl;

        //update the gui
        //resize it manally as it doesn't work atm
        chart->axisX()->setRange(xmin*1.2,xmax*1.2);
        chart->axisY()->setRange(ymin*1.2,ymax*1.2);
        app.processEvents();

        //timer
        usleep(100000);
    }

}
