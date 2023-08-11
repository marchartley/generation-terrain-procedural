#ifndef DISPLAYGRAPHICS_H
#define DISPLAYGRAPHICS_H

#include "DataStructure/Matrix3.h"
#include "DataStructure/Vector3.h"
#include<QtCharts>
#include<QChartView>
#include<QLineSeries>
#include <QPixmap>
#include <QSizePolicy>
#include <iostream>
//#include <QButton>

enum PlotColor {
    WHITE, GRAY, BLACK, RED, GREEN, BLUE, RANDOM
};

inline std::map<PlotColor, QColor> PlotColorToQColor = {
    {WHITE, Qt::white},
    {GRAY, Qt::gray},
    {BLACK, Qt::black},
    {RED, Qt::red},
    {GREEN, Qt::green},
    {BLUE, Qt::blue}
};

class Chart;
class ChartView : public QChartView {
    Q_OBJECT
public:
    ChartView(QWidget *parent = nullptr);
    ChartView(QChart *chart, QWidget *parent = nullptr);
    ChartView(Chart *chart, QWidget *parent = nullptr);

    void lockView() { this->locked = true; }
    void unlockView() { this->locked = false; }

    QPoint previousMousePos;
protected:
    bool viewportEvent(QEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void keyPressEvent(QKeyEvent *event);
    void wheelEvent(QWheelEvent *event);

    bool locked = false;

public:
    Chart* _chart;

Q_SIGNALS:
    void updated();
    void clickedOnValue(Vector3 pos);
};

class Chart : public QChart {
    Q_OBJECT
public:
    Chart();

protected:
    bool sceneEvent(QEvent *event);
    bool gestureEvent(QGestureEvent *event);
//    virtual void wheelEvent(QWheelEvent* event) override;
};

class Plotter : public QDialog {
    Q_OBJECT
private: // Singleton
    Plotter(QWidget* parent = nullptr);
    Plotter(ChartView* chartView, QWidget* parent = nullptr);

public:
    static Plotter* getInstance();
    static Plotter* init(ChartView* chartView = nullptr, QWidget* parent = nullptr);

    void addPlot(std::vector<float> data, std::string name = "", QColor color = Qt::gray);
    void addPlot(std::vector<Vector3> data, std::string name = "", QColor color = Qt::gray);

    void addScatter(std::vector<float> data, std::string name = "", std::vector<std::string> labels = std::vector<std::string>(), std::vector<QColor> colors = std::vector<QColor>());
    void addScatter(std::vector<Vector3> data, std::string name = "", std::vector<std::string> labels = std::vector<std::string>(), std::vector<QColor> colors = std::vector<QColor>());

    void addImage(GridV3 image, bool normalize = false);
    void addImage(GridF image, bool normalize = false);
    void addImage(Matrix3<double> image, bool normalize = false);
    void addImage(GridI image, bool normalize = false);

    void draw();
    int exec();
    void saveFig(std::string filename);
    void resizeEvent(QResizeEvent* event);
    void showEvent(QShowEvent* event);

    void reset();

    QPushButton* saveButton;
    ChartView* chartView;
    QImage* backImage = nullptr;
    std::string title;
    std::vector<std::vector<Vector3>> plot_data;
    std::vector<std::string> plot_names;
    std::vector<QColor> plot_colors;
    std::vector<std::vector<Vector3>> scatter_data;
    std::vector<std::vector<std::string>> scatter_labels;
    std::vector<std::vector<QColor>> scatter_colors;
    std::vector<std::string> scatter_names;

    std::vector<std::vector<QGraphicsTextItem*>> graphicLabels;

    std::vector<std::pair<int, int>> selectedScatterData;
    std::vector<std::pair<int, int>> selectedPlotData;


private:
    static Plotter* instance;
//    QValueAxis* m_axisX;
//    QValueAxis* m_axisY;
public Q_SLOTS:
    void updateLabelsPositions();
    void selectData(Vector3 pos);
};


class TextItem : public QGraphicsItem {
public:
    TextItem(QString text, QPoint pos, QChart *chart, QAbstractSeries *series);
    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) override;

    void setText(const QString &text);
    void setAnchor(QPointF point);
private:
    QChart *_chart;
    QAbstractSeries *_series;
    QString _text;
    QRectF _textRect;
    QPointF _anchor;
};

#endif // DISPLAYGRAPHICS_H
