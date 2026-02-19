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
#include <vector>

#include "GUIElements/CommonInterface.h"
#include "Utils/Utils.h"
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


    ChartView* setOverlay(const GridV3& image, std::string layerName = "default", const GridF& alpha = GridF(1, 1, 1, 1.f), int overlayLayer = 0);
    // ChartView* setOverlay(const GridF& image, std::string layerName = "default", const GridF& alpha = GridF(1, 1, 1, 1.f));
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
    Chart* _chart = nullptr;
    PlotModel* _dataModel = nullptr;

Q_SIGNALS:
    void updated();
    void clickedOnValue(const Vector3& pos, bool leftClick, bool rightClick);
    void mouseMoved(const Vector3& relativePos, const Vector3& previousMousePos, QMouseEvent* e);
    void keyPressed(QKeyEvent* event);
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

    PlotImageData* setImage(const GridV3& img);
    PlotImageData* setImage(const GridF& img);
    PlotImageData* setNormalized(bool normalize);
    PlotImageData* setColorRanges(const Vector3& minRange, const Vector3& maxRange);
    PlotImageData* setAbsolute(bool absolute);
    PlotImageData* setClamped(bool clamp);

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

    PlotVectorData* setField(const GridV3 &field);

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

    PlotModel* addPlot(const std::vector<Vector3>& data, const std::string& name = "", const QColor& color = Qt::gray);

    PlotModel* addScatter(const std::vector<Vector3>& data, const std::string& name = "", const std::vector<std::string>& labels = std::vector<std::string>(), std::vector<QColor> colors = std::vector<QColor>());

    PlotModel* addImage(const GridV3& image);
    PlotModel* addImage(const GridF& image);

    PlotModel* addVectorField(const GridV3& field);

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
    PlotVectorData vectorData;

    GridV3 getImage() const { return imageData.getImage(); }
    GridF getImageGrey() const { return imageData.getImageGrey(); }
    // GridV3 displayedImage;
    // QImage* backImage = nullptr;

    std::string title;
};




#define DECLARE_PLOTTER_GETTER(Type) \
public: \
static Type* get(const std::string& name = "") { return AbstractPlotter::getInstance<Type>([](const std::string& name, ChartView* cv, QWidget* p){ return new Type(name, cv, p); }, toLower(name)); }


class AbstractPlotter : public QDialog {
    Q_OBJECT
protected: // Singleton
    AbstractPlotter(const std::string& name, const std::string& title = "", QWidget* parent = nullptr);
    AbstractPlotter(const std::string& name, ChartView* chartView, const std::string& title = "", QWidget* parent = nullptr);

    template <class Derive>
    static std::string getIDname(const std::string& name) { return name + typeid(Derive).name(); }

public:
    template <class Derived, class Factory>
    static Derived* getInstance(Factory&& factory, std::string name = "") {
        if (name == "") name = Derived::defaultName;
        if (Derived::instances.count(getIDname<Derived>(name)) == 0) {
            init<Derived>(factory, name);
        }
        return dynamic_cast<Derived*>(Derived::instances[getIDname<Derived>(name)]);
    }
    template <class Derived, class Factory>
    static void init(Factory&& factory, const std::string& name, ChartView* chartView = nullptr, QWidget* parent = nullptr) {
        if (Derived::instances.count(getIDname<Derived>(name)))
            delete Derived::instances[getIDname<Derived>(name)];
        Derived::instances[getIDname<Derived>(name)] = factory(name, chartView, parent);
    }
    // static AbstractPlotter* getInstance(std::string name = "");
    // static AbstractPlotter* get(std::string name = "") { return AbstractPlotter::getInstance(toLower(name)); }
    // static AbstractPlotter* init(const std::string& name, ChartView* chartView = nullptr, QWidget* parent = nullptr);

    AbstractPlotter* addPlot(std::vector<float> data, std::string name = "", QColor color = Qt::gray);
    AbstractPlotter* addPlot(std::vector<Vector3> data, std::string name = "", QColor color = Qt::gray);
    AbstractPlotter* addPlot(const BSpline& data, std::string name = "", QColor color = Qt::gray);

    AbstractPlotter* addScatter(std::vector<float> data, std::string name = "", std::vector<std::string> labels = std::vector<std::string>(), std::vector<QColor> colors = std::vector<QColor>());
    AbstractPlotter* addScatter(std::vector<Vector3> data, std::string name = "", std::vector<std::string> labels = std::vector<std::string>(), std::vector<QColor> colors = std::vector<QColor>());

    AbstractPlotter* addImage(const GridV3 &image);
    AbstractPlotter* addImage(const GridF& image);
    AbstractPlotter* addImage(const Matrix3<double>& image);
    AbstractPlotter* addImage(const GridI& image);

    AbstractPlotter* addVectorField(const GridV3& field);

    AbstractPlotter* setOverlay(const GridV3& colors, const GridF& alpha, const std::string& overlayName = "default", int layer = 0);
    AbstractPlotter* setOverlay(const std::pair<GridV3, GridF>& colorAndAlpha, const std::string& overlayName = "default", int layer = 0);
    AbstractPlotter* setOverlay(const GridF& colors, const GridF& alpha, const std::string& overlayName = "default", int layer = 0);
    AbstractPlotter* setOverlay(const std::pair<GridF, GridF>& colorAndAlpha, const std::string& overlayName = "default", int layer = 0);
    AbstractPlotter* showOverlay(const std::string& overlayName = "default");
    AbstractPlotter* hideOverlay(const std::string& overlayName = "default");

    GridV3 computeStreamLinesRendering(const GridV3& field, Vector3 imgSize = Vector3::invalid) const;
    AbstractPlotter* addStreamLines(const GridV3& field, Vector3 imgSize = Vector3::invalid, float opacity = .5f);

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
    bool hasVectorField() const { return !this->dataModel->vectorData.field.empty(); }

    ChartView* chartView;
    PlotModel* dataModel;
    QLabel* mouseInfoLabel = nullptr;

    InterfaceUI* mainInterface;
    InterfaceUI* toolsInterface;
    InterfaceUI* viewOptionsInterface;
    InterfaceUI* saveCopyInterface;
    InterfaceUI* infosInterface;

    InterfaceUI* viewAndCopyInterface;

    std::string name;


protected:
    static std::string defaultName;
    static std::map<std::string, AbstractPlotter*> instances;
    static std::string id_name;
    //    QValueAxis* m_axisX;
    //    QValueAxis* m_axisY;
public Q_SLOTS:
    // AbstractPlotter* updateLabelsPositions();
    // AbstractPlotter* selectData(const Vector3& pos);
    virtual AbstractPlotter* displayInfoUnderMouse(const Vector3& relativeMousePos);
    virtual AbstractPlotter* draw();
    virtual AbstractPlotter* show();
    virtual AbstractPlotter* updateUI(bool forceUpdate = false);


    virtual AbstractPlotter* updateToolsInterface();
    virtual AbstractPlotter* updateViewOptionsInterface();
    virtual AbstractPlotter* updateSaveCopyInterface();
    virtual AbstractPlotter* updateInfosInterface();

Q_SIGNALS:
    void clickedOnImage(const Vector3& pos, Vector3 value, bool leftClick, bool rightClick);
    void movedOnImage(const Vector3& pos, const Vector3& previousPos, QMouseEvent* event);
    void updated();
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









template <class T>
Matrix3<T>& PlotVectorData::drawLine(Matrix3<T>& img, const T& color, const Vector3& start, const Vector3& end) {
    auto line = (end - start);
    int dx = line.x();
    int dy = line.y();

    // calculate steps required for generating pixels
    int steps = abs(dx) > abs(dy) ? abs(dx) : abs(dy);

    // calculate increment in x & y for each steps
    float Xinc = dx / (float)steps;
    float Yinc = dy / (float)steps;

    // Put pixel for each step
    auto p = start;
    for (int i = 0; i <= steps; i++) {
        img[p] = color;
        p.x() += Xinc;
        p.y() += Yinc;
    }

    return img;
}



template<class T>
Matrix3<T> &PlotVectorData::drawCircle(Matrix3<T> &img, const T &color, const Vector3 &center, float radius)
{
    // Function to put pixels
    // at subsequence points
    auto drawCircle = [&](int xc, int yc, int x, int y){
        img(Vector3i(xc+x, yc+y)) = color;
        img(Vector3i(xc-x, yc+y)) = color;
        img(Vector3i(xc+x, yc-y)) = color;
        img(Vector3i(xc-x, yc-y)) = color;
        img(Vector3i(xc+y, yc+x)) = color;
        img(Vector3i(xc-y, yc+x)) = color;
        img(Vector3i(xc+y, yc-x)) = color;
        img(Vector3i(xc-y, yc-x)) = color;
    };

    int xc = center.x();
    int yc = center.y();
    int r = radius;
    // Function for circle-generation
    // using Bresenham's algorithm
    int x = 0, y = r;
    int d = 3 - 2 * r;
    drawCircle(xc, yc, x, y);
    while (y >= x){

        // check for decision parameter
        // and correspondingly
        // update d, y
        if (d > 0) {
            y--;
            d = d + 4 * (x - y) + 10;
        }
        else
            d = d + 4 * x + 6;

        // Increment x after updating decision parameter
        x++;

        // Draw the circle using the new coordinates
        drawCircle(xc, yc, x, y);
    }
    return img;
}

#endif // ABSTRACTPLOTTER_H
