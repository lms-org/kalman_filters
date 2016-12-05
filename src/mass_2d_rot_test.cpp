#include <kalman_filter/ctrv_vxy.h>
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
#include <QtWidgets>
#include <filter/adam.h>

Eigen::Vector2f vV;
Eigen::Vector2f xV;

void scaleChart(QtCharts::QChart *chart, std::vector<QtCharts::QXYSeries*> series_,const bool sameAxes){
    float xMin = std::numeric_limits<float>::max(); // everything is <= this
    float xMax = std::numeric_limits<float>::min(); // everything is >= this
    float yMin = std::numeric_limits<float>::max();
    float yMax = std::numeric_limits<float>::min();
    for(QtCharts::QXYSeries* series:series_){
        for (QPointF &p: series->points()) {
            xMin = qMin<float>(xMin, p.x());
            xMax = qMax<float>(xMax, p.x());
            yMin = qMin<float>(yMin, p.y());
            yMax = qMax<float>(yMax, p.y());
        }
    }
    if(sameAxes){
        xMin = qMin<float>(xMin, yMin);
        xMax = qMin<float>(xMax, yMax);
        yMin = xMin;
        yMax = xMax;
    }
    chart->axisX()->setRange(xMin-1,xMax+1);
    chart->axisY()->setRange(yMin-1,yMax+1);

}

void scaleChart(QtCharts::QChart *chart, QtCharts::QXYSeries *series, const bool sameAxes){
    float xMin = std::numeric_limits<float>::max(); // everything is <= this
    float xMax = std::numeric_limits<float>::min(); // everything is >= this
    float yMin = std::numeric_limits<float>::max();
    float yMax = std::numeric_limits<float>::min();
    for (QPointF &p: series->points()) {
        xMin = qMin<float>(xMin, p.x());
        xMax = qMax<float>(xMax, p.x());
        yMin = qMin<float>(yMin, p.y());
        yMax = qMax<float>(yMax, p.y());
    }
    if(sameAxes){
        xMin = qMin<float>(xMin, yMin);
        xMax = qMax<float>(xMax, yMax);
        yMin = xMin;
        yMax = xMax;
    }
    chart->axisX()->setRange(xMin-1,xMax+1);
    chart->axisY()->setRange(yMin-1,yMax+1);
}

int main(int argc, char** argv){
    vV(0) = 1;
    vV(1) = 0;
    xV.setZero();
    kalman_filters::ctrv_vxy::MassModelUKF ukf;
    ukf.init();
    //QWT http://qwt.sourceforge.net/
    //oder einfach nur qt5 http://doc.qt.io/qt-5/qtdatavisualization-index.html

    QT_CHARTS_USE_NAMESPACE

    QApplication app(argc, argv);
    QWidget *central = new QWidget();
    QGridLayout *layout = new QGridLayout;
    central->setLayout(layout);

    QLineSeries ukfVY,ukfVX,ukfV;
    QLineSeries points;

    //velocity X
    QChart *vxChart = new QChart();
    ukfVY.setName("ukf vy");
    ukfVX.setName("ukf vx");
    ukfV.setName("ukf abs v");
    vxChart->addSeries(&ukfVY);
    vxChart->addSeries(&ukfVX);
    vxChart->addSeries(&ukfV);
    vxChart->createDefaultAxes();
    vxChart->setTitle("UKF");
    QChartView vxChartView(vxChart);
    vxChartView.setRenderHint(QPainter::Antialiasing);
    layout->addWidget(&vxChartView);


    QChart *positionChart = new QChart();
    points.setName("position");
    positionChart->addSeries(&points);
    positionChart->createDefaultAxes();
    positionChart->legend()->hide();
    vxChart->setTitle("Position in m");
    //velocity X
    QChartView positionChartView(positionChart);
    positionChartView.setRenderHint(QPainter::Antialiasing);
    layout->addWidget(&positionChartView);

    QMainWindow window;
    window.setCentralWidget(central);

    window.resize(400, 300);
    window.show();
    float velocityX = 5;
    float velocityY = 0;
    float angularVelocity = 0.2;
    int iter = 0;

    while(true){
        //std::cout<<"vx: "<<velocityX<<std::endl;
        ukf.setMeasurementVec(velocityX,velocityY,angularVelocity);
        //ukf.predict(0.1); //TODO we just test predict
        ukf.pu(0.1);
        /*
        Eigen::Matrix2f vRot;
        vRot(0,0) = std::cos(0.1);
        vRot(0,1) = -std::sin(0.1);
        vRot(1,0) = std::sin(0.1);
        vRot(1,1) = std::cos(0.1);
        vV = vRot*vV;
        xV += vV*0.1;

        ukfV.append(iter,std::sqrt(vV(0)*vV(0)+vV(1)*vV(1)));
        ukfVX.append(iter,vV(0));
        ukfVY.append(iter,vV(1));//ukf.lastState.vx());
        points.append(xV(0),xV(1));
        */

        ukfVX.append(iter,ukf.lastState.vx());
        ukfVY.append(iter,ukf.lastState.vy());
        ukfV.append(iter,std::sqrt(ukf.lastState.vx()*ukf.lastState.vx()+ukf.lastState.vy()*ukf.lastState.vy()));
        points.append(ukf.lastState.x(),ukf.lastState.y());
        //std::cout<<"vx ukf: "<<ukf.lastState.vx()<<std::endl;

        iter++;

        //update the gui
        //resize it manally as it doesn't work atm
        std::vector<QtCharts::QXYSeries*> series;
        series.push_back(&ukfVY);
        series.push_back(&ukfVX);
        scaleChart(vxChart,series,false);
        scaleChart(positionChart,&points,true);
        app.processEvents();

        //timer
        usleep(100000);
    }

}
