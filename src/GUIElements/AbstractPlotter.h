#ifndef ABSTRACTPLOTTER_H
#define ABSTRACTPLOTTER_H

#include "DataStructure/Matrix3.h"
#include "DataStructure/Vector3.h"
// #include<QtCharts>
// #include<QChartView>
// #include<QLineSeries>
// #include <QPixmap>
// #include <QSizePolicy>
// #include <iostream>
//#include <QButton>
#include <vector>

#include "GUIElements/CommonInterface.h"
#include "Utils/Utils.h"
#include "DataStructure/Image.h"
#include "ChartView.h"


#define DECLARE_PLOTTER_GETTER(Type) \
public: \
static Type* get(const std::string& name = "") { return AbstractPlotter::getInstance<Type>([](const std::string& name, ChartView* cv, QWidget* p){ return new Type(name, cv, p); }, toLower(name)); } \
static Type* reset(const std::string& name = "") { AbstractPlotter::init<Type>([](const std::string& name, ChartView* cv, QWidget* p){ return new Type(name, cv, p); }, toLower(name)); return Type::get(name);}


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
            // delete Derived::instances[getIDname<Derived>(name)];
            Derived::instances[getIDname<Derived>(name)]->deleteLater();
        // AbstractPlotter::erase<Derived>(getIDname<Derived>(name));
        Derived::instances.erase(getIDname<Derived>(name));
        Derived::instances[getIDname<Derived>(name)] = factory(name, chartView, parent);
    }

    template <class Derived>
    static void erase(const std::string& name) {
        if (Derived::instances.count(getIDname<Derived>(name)) > 0) delete Derived::instances[getIDname<Derived>(name)];
    }
    // static AbstractPlotter* getInstance(const std::string& name = "");
    // static AbstractPlotter* get(const std::string& name = "") { return AbstractPlotter::getInstance(toLower(name)); }
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
    AbstractPlotter* saveFig(const std::string& filename);
    AbstractPlotter* copyToClipboard();
    void resizeEvent(QResizeEvent* event);
    void showEvent(QShowEvent* event);
    void hideEvent(QHideEvent *event);

    QTimer *animate(std::function<void()> callback, int interval_ms = 30);

    // AbstractPlotter* reset();

    bool hasPlotValues() const { return !this->dataModel->plot_data.empty(); }
    bool hasScatterValues() const { return !this->dataModel->scatter_data.empty(); }
    bool hasImage() const { return !this->dataModel->getImage().empty(); }
    bool hasVectorField() const { return !this->dataModel->vectorData.field.empty(); }

    ChartView* chartView = nullptr;
    PlotModel* dataModel = nullptr;
    QLabel* mouseInfoLabel = nullptr;

    InterfaceUI* mainInterface = nullptr;
    InterfaceUI* toolsInterface = nullptr;
    InterfaceUI* viewOptionsInterface = nullptr;
    InterfaceUI* saveCopyInterface = nullptr;
    InterfaceUI* infosInterface = nullptr;

    InterfaceUI* viewAndCopyInterface = nullptr;

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
