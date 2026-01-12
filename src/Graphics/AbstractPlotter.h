#ifndef ABSTRACTPLOTTER_H
#define ABSTRACTPLOTTER_H

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


    ChartView* setOverlay(GridV3 image, GridF alpha = GridF(1, 1, 1, 1.f));
    GridV3 overlayColors;
    GridF  overlayAlpha;

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

class PlotImageData {
public:
    PlotImageData();
    PlotImageData(const GridV3& img);

    PlotImageData* setImage(const GridV3& img);
    PlotImageData* setNormalized(bool normalize);
    PlotImageData* setColorRanges(const Vector3& minRange, const Vector3& maxRange);
    PlotImageData* setAbsolute(bool absolute);
    PlotImageData* setClamped(bool clamp);

    GridV3& getImage() { return this->image; }
    const GridV3& getImage() const { return this->image; }
    QImage computeDisplayedImage() const;
    QImage computeDisplayedImage(const GridV3& overlay, const GridF& overlayAlpha) const;


// protected:
    bool normalized = false;
    bool absolute = false;
    bool clamped = true;
    Vector3 colorRangeMin = Vector3::min();
    Vector3 colorRangeMax = Vector3::max();
    Vector3i displayedColors = Vector3i(1, 1, 1);
    GridV3 image;
};

class PlotModel {
public:
    PlotModel();

    PlotModel* addPlot(const std::vector<Vector3>& data, const std::string& name = "", const QColor& color = Qt::gray);

    PlotModel* addScatter(const std::vector<Vector3>& data, const std::string& name = "", const std::vector<std::string>& labels = std::vector<std::string>(), std::vector<QColor> colors = std::vector<QColor>());

    PlotModel* addImage(const GridV3& image/*, bool clamped = false, bool normalized = false, bool absolute = false, const Vector3& minColors = Vector3::min(), const Vector3& maxColors = Vector3::max()*/);

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

    PlotImageData imageData;

    GridV3& getImage() { return imageData.getImage(); }
    const GridV3& getImage() const { return imageData.getImage(); }
    // GridV3 displayedImage;
    // QImage* backImage = nullptr;

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

    AbstractPlotter* setOverlay(const GridV3& colors, const GridF& alpha);

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

    bool hasPlotValues() const { return !this->dataModel->plot_data.empty(); }
    bool hasScatterValues() const { return !this->dataModel->scatter_data.empty(); }
    bool hasImage() const { return !this->dataModel->getImage().empty(); }

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

#endif // ABSTRACTPLOTTER_H
