#ifndef PLOTTER_H
#define PLOTTER_H

#include "DataStructure/Matrix3.h"
#include "DataStructure/Vector3.h"
#include<QtCharts>
#include<QChartView>
#include<QLineSeries>
#include <QPixmap>
#include <QSizePolicy>
#include <iostream>
//#include <QButton>

#include "Interface/CommonInterface.h"
#include "Utils/Utils.h"

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
class PlotModel;
class ChartView : public QChartView {
    Q_OBJECT
public:
    ChartView(QWidget *parent = nullptr);
    ChartView(QChart *chart, QWidget *parent = nullptr);
    ChartView(Chart *chart, QWidget *parent = nullptr);

    void lockView() { this->locked = true; }
    void unlockView() { this->locked = false; }

    ChartView* setPlotModel(PlotModel* dataModel, std::string title = "");
    ChartView* updateLabelsPositions();

    bool selectData(const Vector3& pos);

    Vector3 getRelativeMousePositionInImage(const Vector3& pos);

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
    Chart* _chart = nullptr;
    PlotModel* _dataModel = nullptr;

Q_SIGNALS:
    void updated();
    void clickedOnValue(const Vector3& pos);
    void mouseMoved(const Vector3& relativePos, const Vector3& previousMousePos, QMouseEvent* e);
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

class PlotModel {
public:
    PlotModel();

    PlotModel* addPlot(std::vector<Vector3> data, std::string name = "", QColor color = Qt::gray);

    PlotModel* addScatter(std::vector<Vector3> data, std::string name = "", std::vector<std::string> labels = std::vector<std::string>(), std::vector<QColor> colors = std::vector<QColor>());

    PlotModel* addImage(GridV3 image, bool clamped = false, bool normalized = false, bool absolute = false, Vector3 minColors = Vector3::min(), Vector3 maxColors = Vector3::max());

    PlotModel* reset();

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

    GridV3 displayedImage;
    QImage* backImage = nullptr;

    std::string title;
};







class AbstractPlotter : public QDialog {
    Q_OBJECT
protected: // Singleton
    AbstractPlotter(std::string name, QWidget* parent = nullptr);
    AbstractPlotter(std::string name, ChartView* chartView, QWidget* parent = nullptr);

public:
    // static AbstractPlotter* getInstance(std::string name = "");
    // static AbstractPlotter* get(std::string name = "") { return AbstractPlotter::getInstance(toLower(name)); }
    // static AbstractPlotter* init(std::string name, ChartView* chartView = nullptr, QWidget* parent = nullptr);

    AbstractPlotter* addPlot(std::vector<float> data, std::string name = "", QColor color = Qt::gray);
    AbstractPlotter* addPlot(std::vector<Vector3> data, std::string name = "", QColor color = Qt::gray);
    AbstractPlotter* addPlot(const BSpline& data, std::string name = "", QColor color = Qt::gray);

    AbstractPlotter* addScatter(std::vector<float> data, std::string name = "", std::vector<std::string> labels = std::vector<std::string>(), std::vector<QColor> colors = std::vector<QColor>());
    AbstractPlotter* addScatter(std::vector<Vector3> data, std::string name = "", std::vector<std::string> labels = std::vector<std::string>(), std::vector<QColor> colors = std::vector<QColor>());

    AbstractPlotter* addImage(GridV3 image);
    AbstractPlotter* addImage(const GridF& image);
    AbstractPlotter* addImage(const Matrix3<double>& image);
    AbstractPlotter* addImage(const GridI& image);

    GridV3 computeVectorFieldRendering(const GridV3& field, float reductionFactor = .1f, Vector3 imgSize = Vector3(false)) const;
    AbstractPlotter* addVectorField(const GridV3& field, float reductionFactor = .1f, Vector3 imgSize = Vector3(false), float opacity = .5f);
    GridV3 computeStreamLinesRendering(const GridV3& field, Vector3 imgSize = Vector3(false)) const;
    AbstractPlotter* addStreamLines(const GridV3& field, Vector3 imgSize = Vector3(false), float opacity = .5f);

    int exec();
    AbstractPlotter* saveFig(std::string filename);
    AbstractPlotter* copyToClipboard();
    void resizeEvent(QResizeEvent* event);
    void showEvent(QShowEvent* event);

    QTimer *animate(std::function<void()> callback, int interval_ms = 30);

    AbstractPlotter* reset();

    ChartView* chartView;
    PlotModel* dataModel;
    QLabel* mouseInfoLabel = nullptr;

    InterfaceUI* toolsInterface;
    InterfaceUI* viewOptionsInterface;
    InterfaceUI* saveCopyInterface;
    InterfaceUI* infosInterface;

    std::string name;

protected:
    static std::string defaultName;
    static std::map<std::string, AbstractPlotter*> instances;
    //    QValueAxis* m_axisX;
    //    QValueAxis* m_axisY;
public Q_SLOTS:
    // AbstractPlotter* updateLabelsPositions();
    // AbstractPlotter* selectData(const Vector3& pos);
    AbstractPlotter* displayInfoUnderMouse(const Vector3& relativeMousePos);
    AbstractPlotter* draw();
    AbstractPlotter* show();
    AbstractPlotter* updateUI();

Q_SIGNALS:
    void clickedOnImage(const Vector3& pos, Vector3 value);
    void movedOnImage(const Vector3& pos, const Vector3& previousPos, QMouseEvent* event);
};





class Plotter : public AbstractPlotter {
    Q_OBJECT
private: // Singleton
    Plotter(std::string name, QWidget* parent = nullptr);
    Plotter(std::string name, ChartView* chartView, QWidget* parent = nullptr);

public:
    static Plotter* getInstance(std::string name = "");
    static Plotter* get(std::string name = "") { return Plotter::getInstance(toLower(name)); }
    static Plotter* init(std::string name, ChartView* chartView = nullptr, QWidget* parent = nullptr);

    // Plotter* addPlot(std::vector<float> data, std::string name = "", QColor color = Qt::gray);
    // Plotter* addPlot(std::vector<Vector3> data, std::string name = "", QColor color = Qt::gray);
    // Plotter* addPlot(const BSpline& data, std::string name = "", QColor color = Qt::gray);

    // Plotter* addScatter(std::vector<float> data, std::string name = "", std::vector<std::string> labels = std::vector<std::string>(), std::vector<QColor> colors = std::vector<QColor>());
    // Plotter* addScatter(std::vector<Vector3> data, std::string name = "", std::vector<std::string> labels = std::vector<std::string>(), std::vector<QColor> colors = std::vector<QColor>());

    // Plotter* addImage(GridV3 image);
    // Plotter* addImage(const GridF& image);
    // Plotter* addImage(const Matrix3<double>& image);
    // Plotter* addImage(const GridI& image);

    // GridV3 computeVectorFieldRendering(const GridV3& field, float reductionFactor = .1f, Vector3 imgSize = Vector3(false)) const;
    // Plotter* addVectorField(const GridV3& field, float reductionFactor = .1f, Vector3 imgSize = Vector3(false), float opacity = .5f);
    // GridV3 computeStreamLinesRendering(const GridV3& field, Vector3 imgSize = Vector3(false)) const;
    // Plotter* addStreamLines(const GridV3& field, Vector3 imgSize = Vector3(false), float opacity = .5f);

    // int exec();
    // Plotter* saveFig(std::string filename);
    // Plotter* copyToClipboard();
    // void resizeEvent(QResizeEvent* event);
    // void showEvent(QShowEvent* event);

    // QTimer *animate(std::function<void()> callback, int interval_ms = 30);

    // Plotter* reset();

    bool normalizedMode = false;
    bool absoluteMode = false;
    float minValueToDisplay = -10000.f;
    float maxValueToDisplay = 10000.f;
    bool clampValues = false;

    bool displayR = true;
    bool displayG = true;
    bool displayB = true;

    RangeSliderElement* rangeValuesWidget;

public Q_SLOTS:
    Plotter* updateUI();

    Plotter* setNormalizedModeImage(bool normalize);
    Plotter* setAbsoluteModeImage(bool absolute);
    Plotter* setFilteredValuesImage(bool filtered);

// Q_SIGNALS:
    // void clickedOnImage(const Vector3& pos, Vector3 value);
    // void movedOnImage(const Vector3& pos, const Vector3& previousPos, QMouseEvent* event);
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









/*

class Plotter : public QDialog {
    Q_OBJECT
private: // Singleton
    Plotter(std::string name, QWidget* parent = nullptr);
    Plotter(std::string name, ChartView* chartView, QWidget* parent = nullptr);

public:
    static Plotter* getInstance(std::string name = "");
    static Plotter* get(std::string name = "") { return Plotter::getInstance(toLower(name)); }
    static Plotter* init(std::string name, ChartView* chartView = nullptr, QWidget* parent = nullptr);

    Plotter* addPlot(std::vector<float> data, std::string name = "", QColor color = Qt::gray);
    Plotter* addPlot(std::vector<Vector3> data, std::string name = "", QColor color = Qt::gray);
    Plotter* addPlot(const BSpline& data, std::string name = "", QColor color = Qt::gray);

    Plotter* addScatter(std::vector<float> data, std::string name = "", std::vector<std::string> labels = std::vector<std::string>(), std::vector<QColor> colors = std::vector<QColor>());
    Plotter* addScatter(std::vector<Vector3> data, std::string name = "", std::vector<std::string> labels = std::vector<std::string>(), std::vector<QColor> colors = std::vector<QColor>());

    Plotter* addImage(GridV3 image);
    Plotter* addImage(const GridF& image);
    Plotter* addImage(const Matrix3<double>& image);
    Plotter* addImage(const GridI& image);

    GridV3 computeVectorFieldRendering(const GridV3& field, float reductionFactor = .1f, Vector3 imgSize = Vector3(false)) const;
    Plotter* addVectorField(const GridV3& field, float reductionFactor = .1f, Vector3 imgSize = Vector3(false), float opacity = .5f);
    GridV3 computeStreamLinesRendering(const GridV3& field, Vector3 imgSize = Vector3(false)) const;
    Plotter* addStreamLines(const GridV3& field, Vector3 imgSize = Vector3(false), float opacity = .5f);

    int exec();
    Plotter* saveFig(std::string filename);
    Plotter* copyToClipboard();
    void resizeEvent(QResizeEvent* event);
    void showEvent(QShowEvent* event);

    QTimer *animate(std::function<void()> callback, int interval_ms = 30);

    Plotter* reset();

    bool normalizedMode = false;
    bool absoluteMode = false;
    float minValueToDisplay = -10000.f;
    float maxValueToDisplay = 10000.f;
    bool clampValues = false;

    bool displayR = true;
    bool displayG = true;
    bool displayB = true;

    RangeSliderElement* rangeValuesWidget;
//    ButtonElement* saveButton;
    ChartView* chartView;
    QImage* backImage = nullptr;
    QLabel* mouseInfoLabel = nullptr;
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

    GridV3 displayedImage;

    InterfaceUI* interfaceButtons;

    std::string name;

private:
    static std::string defaultName;
    static std::map<std::string, Plotter*> instances;
//    QValueAxis* m_axisX;
//    QValueAxis* m_axisY;
public Q_SLOTS:
    Plotter* updateLabelsPositions();
    Plotter* selectData(const Vector3& pos);
    Plotter* displayInfoUnderMouse(const Vector3& relativeMousePos);
    Plotter* draw();
    Plotter* show();
    Plotter* updateUI();

    Plotter* setNormalizedModeImage(bool normalize);
    Plotter* setAbsoluteModeImage(bool absolute);
    Plotter* setFilteredValuesImage(bool filtered);

Q_SIGNALS:
    void clickedOnImage(const Vector3& pos, Vector3 value);
    void movedOnImage(const Vector3& pos, const Vector3& previousPos, QMouseEvent* event);
};
*/

#endif // PLOTTER_H
