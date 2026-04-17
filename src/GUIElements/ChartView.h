#ifndef CHARTVIEW_H
#define CHARTVIEW_H

#include <QColor>
#include <QtCharts>
#include <QChartView>

#include "GUIElements/PlottingData.h"
#include "GUIElements/ImageData.h"
#include "GUIElements/VectorFieldData.h"

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

class PlotModel;

class Chart : public QChart {
    Q_OBJECT
public:
    Chart(QGraphicsItem* parent = nullptr);

protected:
    bool sceneEvent(QEvent *event);
    bool gestureEvent(QGestureEvent *event);
    //    virtual void wheelEvent(QWheelEvent* event) override;
};

class ChartView : public QChartView {
    Q_OBJECT
public:
    ChartView(QWidget *parent = nullptr);

    void lockView() { this->locked = true; }
    void unlockView() { this->locked = false; }

    ChartView& setPlotModel(std::shared_ptr<PlotModel> dataModel, const std::string& title = "");
    ChartView& updateLabelsPositions();

    bool selectData(const Vector3& pos);

    Vector3 getRelativeMousePositionInImage(const Vector3& pos);

    QPoint previousMousePos;


    const ChartView& setOverlay(const GridV3& image, std::string layerName = "default", const GridF& alpha = GridF(1, 1, 1, 1.f), int overlayLayer = 0);
    // const ChartView& setOverlay(const GridF& image, std::string layerName = "default", const GridF& alpha = GridF(1, 1, 1, 1.f));
    std::map<std::string, GridV3> overlayColors;
    std::map<std::string, GridF>  overlayAlpha;
    std::map<std::string, bool>  overlayDisplayed;
    std::map<std::string, int>  overlayLayer;

protected:
    bool viewportEvent(QEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void keyPressEvent(QKeyEvent *event);
    void wheelEvent(QWheelEvent *event);

    bool locked = false;

public:
    std::shared_ptr<PlotModel> _dataModel;

Q_SIGNALS:
    void updated();
    void clickedOnValue(const Vector3& pos, bool leftClick, bool rightClick);
    void mouseMoved(const Vector3& relativePos, const Vector3& previousMousePos, QMouseEvent* e);
    void keyPressed(QKeyEvent* event);
};




class PlotModel {
public:
    PlotModel();

    PlotModel& addPlot(const std::vector<Vector3>& data, const std::string& name = "", const QColor& color = Qt::gray);

    PlotModel& addScatter(const std::vector<Vector3>& data, const std::string& name = "", const std::vector<std::string>& labels = std::vector<std::string>(), std::vector<QColor> colors = std::vector<QColor>());

    PlotModel& addImage(const GridV3& image);
    PlotModel& addImage(const GridF& image);

    PlotModel& addVectorField(const GridV3& field);

    PlotModel& reset();

    std::vector<std::pair<int, int>> selectedScatterData;
    std::vector<std::pair<int, int>> selectedPlotData;

    PlotImageData imageData;
    PlotVectorData vectorData;
    PlotLineData plotLineData;
    PlotScatterData scatterData;

    GridV3 getImage() const { return imageData.getImage(); }
    GridF getImageGrey() const { return imageData.getImageGrey(); }

    std::string title;
};

#endif // CHARTVIEW_H
