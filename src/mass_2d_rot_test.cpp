#include <kalman_filter/mass_2d_rot.h>
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
    KALMAN_FILTERS::MASS_2D_ROT_EKF::MassModelUKF ukf;
    ukf.init();
    //QWT http://qwt.sourceforge.net/
    //oder einfach nur qt5 http://doc.qt.io/qt-5/qtdatavisualization-index.html

    QT_CHARTS_USE_NAMESPACE

    QApplication app(argc, argv);

    QLineSeries vx,ukfVX;
    QScatterSeries points;

    QChart *chart = new QChart();
    vx.setName("vx");
    ukfVX.setName("ukf vx");
    chart->addSeries(&vx);
    chart->addSeries(&ukfVX);
    chart->createDefaultAxes();
    chart->setTitle("UKF");

    QChartView chartView(chart);
    chartView.setRenderHint(QPainter::Antialiasing);

    QMainWindow window;
    window.setCentralWidget(&chartView);
    window.resize(400, 300);
    window.show();
    float velocityX = 0;
    float velocityY = 0;
    int iter = 0;
    while(true){
        std::cout<<"vx: "<<velocityX<<std::endl;
        ukf.setMeasurementVec(velocityX,0,0);
        ukf.pu(0.1);
        vx.append(iter,velocityX);
        ukfVX.append(iter,ukf.lastState.vx());
        std::cout<<"vx ukf: "<<ukf.lastState.vx()<<std::endl;


        velocityX += 0.1;
        velocityY += 0.1;
        iter++;

        //std::cout <<"points: "<<mData<<std::endl;

        //update the gui
        //resize it manally as it doesn't work atm
        scaleChart(chart,&vx,false);
        app.processEvents();

        //timer
        usleep(100000);
    }

}
