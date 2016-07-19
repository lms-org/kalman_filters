#include <kalman_filter/mass_ekf.h>
#include <filter/line_line_x.h>
#include <iostream>
#include <Eigen/Eigen>
#include <cmath>


#include <QtGui/QGuiApplication>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtCharts/QChartView>
#include <QtCharts/QSplineSeries>
#include <QtCharts/QScatterSeries>

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


    //QWT http://qwt.sourceforge.net/
    //oder einfach nur qt5 http://doc.qt.io/qt-5/qtdatavisualization-index.html

    QT_CHARTS_USE_NAMESPACE

    QApplication app(argc, argv);

    QLineSeries fit;
    QScatterSeries points;

    QChart *chart = new QChart();
    chart->legend()->hide();
    chart->addSeries(&fit);
    chart->addSeries(&points);
    chart->createDefaultAxes();
    chart->setTitle("Simple line chart example");

    QChartView chartView(chart);
    chartView.setRenderHint(QPainter::Antialiasing);

    QMainWindow window;
    window.setCentralWidget(&chartView);
    window.resize(400, 300);
    window.show();
    //chart->axisX()->setRange(0,10);
    //chart->axisY()->setRange(-1,-1);

    while(true){
        //create new measurement-points
        Eigen::Matrix<double,10,2> mData;
        double d = 1.5;
        double offset = 0;
        float randomF = 1;

        //range
        float xmin = 0;
        float xmax = 0;
        float ymin = 0;
        float ymax = 0;
        for(int i = 0; i < 10; i++){
            float r = ((double) rand() / (RAND_MAX));
            mData(i,0) = i;
            mData(i,1) = d*i+offset + randomF*r;

            //range checks
            if(mData(i,0) < xmin){
                xmin = mData(i,0);
            }
            if(mData(i,0) > xmax){
                xmax = mData(i,0);
            }
            if(mData(i,0) < ymin){
                ymin = mData(i,0);
            }
            if(mData(i,0) > ymax){
                ymax = mData(i,0);
            }
        }
        //conversion to qt points
        points.clear();
        fit.clear();
        for(int i = 0; i < 10; i++){
            points.append(mData(i,0),mData(i,1));
        }

        //TODO fit line

        //update the gui
        //resize it manally as it doesn't work atm
        chart->axisX()->setRange(xmin,xmax);
        chart->axisY()->setRange(ymin,ymax);
        app.processEvents();

        //timer
        usleep(200);
    }

}
