#ifndef CHARTVIEW_H
#define CHARTVIEW_H

#include <QColor>
#include <QtCharts>
#include <QChartView>

#include "DataStructure/Vector3.h"
#include "DataStructure/Matrix3.h"
#include "DataStructure/Image.h"

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




struct DisplayedImageParameters {
    bool normalized = false;
    bool absolute = false;
    bool clamped = true;
    Vector3 colorRangeMin = Vector3::min();
    Vector3 colorRangeMax = Vector3::max();
    Vector3i displayedColors = Vector3i(1, 1, 1);

    BSpline colorRamp = BSpline({Vector3(1, 0, 0), Vector3(1, 1, 1), Vector3(0, 1, 0)});
};

struct PlotImageData {
    // public:
    PlotImageData();
    PlotImageData(const GridV3& img);
    PlotImageData(const GridF& img);

    PlotImageData& setImage(const GridV3& img);
    PlotImageData& setImage(const GridF& img);
    PlotImageData& setNormalized(bool normalize);
    PlotImageData& setColorRanges(const Vector3& minRange, const Vector3& maxRange);
    PlotImageData& setAbsolute(bool absolute);
    PlotImageData& setClamped(bool clamp);

    GridV3 getImage() const { return this->image.getColorImage(); }
    GridF getImageGrey() const { return this->image.getBwImage(); }
    GridV3 prepareImageForDisplay(const Image &img) const;
    QImage computeDisplayedImage(const Vector3i &imgSize = Vector3i::invalid) const;
    QImage computeDisplayedImage(const GridV3& overlay, const GridF& overlayAlpha) const;
    QImage computeDisplayedImage(const std::map<std::string, GridV3>& overlays,
                                 const std::map<std::string, GridF>& overlayAlphas,
                                 const std::map<std::string, bool>& displayedOverlays,
                                 const std::map<std::string, int>& overlayLayers,
                                 const Vector3i& imgSize) const;


    // protected:
    DisplayedImageParameters displayParameters;
    // GridV3 image;
    Image image;
};

struct DisplayedVectorFieldParameters {
    enum VECTOR_DISPLAY { ARROWS, FLOWLINES, NONE };

    VECTOR_DISPLAY displayMode = ARROWS;
    Vector3 backgroundColor = Vector3::white;
    BSpline colorRamp = BSpline({Vector3(70.f, 0.f, 100.f) / 255.f, Vector3(30.f, 160.f, 130.f) / 255.f, Vector3(255.f, 250.f, 0.f)/255.f});
};

struct PlotVectorData {
    PlotVectorData();
    PlotVectorData(const GridV3& field);

    PlotVectorData& setField(const GridV3 &field);

    const GridV3& getField() const { return this->field; }

    std::pair<GridV3, GridF> getFieldImageAndAlpha(const Vector3i &imgSize, const Vector3i &numberOfCells) const;
    static std::pair<GridV3, GridF> createFieldImageAndAlpha(const GridV3& field, Vector3i imgSize, const Vector3i &numberOfCells, DisplayedVectorFieldParameters displayParameters = DisplayedVectorFieldParameters());

    template <class T>
    static Matrix3<T>& drawLine(Matrix3<T>& img, const T& color, const Vector3& start, const Vector3& end);
    template <class T>
    static Matrix3<T>& drawCircle(Matrix3<T>& img, const T& color, const Vector3& center, float radius);

    GridV3 field;

    DisplayedVectorFieldParameters displayParameters;
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
    PlotVectorData vectorData;

    GridV3 getImage() const { return imageData.getImage(); }
    GridF getImageGrey() const { return imageData.getImageGrey(); }
    // GridV3 displayedImage;
    // QImage* backImage = nullptr;

    std::string title;
};

#endif // CHARTVIEW_H
