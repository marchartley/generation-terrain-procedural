#ifndef ABSTRACTPLOTTER_H
#define ABSTRACTPLOTTER_H

#include "DataStructure/Matrix3.h"
#include "DataStructure/Vector3.h"
#include <vector>

#include "GUIElements/CommonInterface.h"
#include "Utils/Utils.h"
#include "DataStructure/Image.h"
#include "GUIElements/ChartView.h"

#include "Utils/Signals.h"

#include "GUIElements/PlottingUtils.h"

#define DECLARE_PLOTTER_GETTER(Type)                    \
public:                                                 \
    static Type& get(const std::string& name = "") {    \
        return AbstractPlotter::getInstance<Type>(      \
           [](const std::string& name, QWidget* p) {    \
                   return new Type(name, p);            \
           },                                           \
           toLower(name));                              \
}                                                       \
    static Type& reset(const std::string& name = "") {  \
        AbstractPlotter::init<Type>(                    \
           [](const std::string& name, QWidget* p) {    \
                   return new Type(name, p);            \
           },                                           \
           toLower(name));                              \
        return Type::get(name);                         \
}

class AbstractPlotter : public QDialog {
    Q_OBJECT
public: // protected: // Singleton
    AbstractPlotter(const std::string& name, QWidget* parent = nullptr);
    AbstractPlotter(const std::string& name, const std::string& title, QWidget* parent = nullptr);
    virtual ~AbstractPlotter() {
    }

    template <class Derive>
    static std::string getIDname(const std::string& name) { return name + typeid(Derive).name(); }

public:
    template <class Derived, class Factory>
    static Derived& getInstance(Factory&& factory, const std::string& _name = "") {
        std::string name = _name;
        if (name == "") name = Derived::defaultName;
        if (Derived::instances.count(getIDname<Derived>(name)) == 0) {
            init<Derived>(factory, name);
        }
        return dynamic_cast<Derived&>(*Derived::instances[getIDname<Derived>(name)]);
    }
    template <class Derived, class Factory>
    static void init(Factory&& factory, const std::string& name, QWidget* parent = nullptr) {
        if (Derived::instances.count(getIDname<Derived>(name)))
            // delete Derived::instances[getIDname<Derived>(name)];
            Derived::instances[getIDname<Derived>(name)]->deleteLater();
        // AbstractPlotter::erase<Derived>(getIDname<Derived>(name));
        Derived::instances.erase(getIDname<Derived>(name));
        Derived::instances[getIDname<Derived>(name)] = factory(name, parent);
    }

    template <class Derived>
    static void erase(const std::string& name) {
        if (Derived::instances.count(getIDname<Derived>(name)) > 0) delete Derived::instances[getIDname<Derived>(name)];
    }
    DECLARE_PLOTTER_GETTER(AbstractPlotter)

    AbstractPlotter& addPlot(const std::vector<float>& data, std::string name = "", Vector3 color = Vector3::gray);
    AbstractPlotter& addPlot(const std::vector<Vector3>& data, std::string name = "", Vector3 color = Vector3::gray);
    AbstractPlotter& addPlot(const Curve& data, std::string name = "", Vector3 color = Vector3::gray);

    AbstractPlotter& addScatter(const std::vector<float>& data, std::string name = "", std::vector<std::string> labels = std::vector<std::string>(), const Vector3& color = Vector3::gray);
    AbstractPlotter& addScatter(const std::vector<float>& dataX, const std::vector<float>& dataY, const std::string& name = "", const std::vector<std::string>& labels = std::vector<std::string>(), const Vector3& color = Vector3::gray);
    AbstractPlotter& addScatter(const std::vector<Vector3>& data, std::string name = "", std::vector<std::string> labels = std::vector<std::string>(), const Vector3& color = Vector3::gray);

    AbstractPlotter& addImage(const GridV3 &image);
    AbstractPlotter& addImage(const GridF& image);
    AbstractPlotter& addImage(const Matrix3<double>& image);
    AbstractPlotter& addImage(const GridI& image);

    AbstractPlotter& addVectorField(const GridV3& field);

    AbstractPlotter& setOverlay(const GridV3& colors, const GridF& alpha, const std::string& overlayName = "default", int layer = 0);
    AbstractPlotter& setOverlay(const std::pair<GridV3, GridF>& colorAndAlpha, const std::string& overlayName = "default", int layer = 0);
    AbstractPlotter& setOverlay(const GridF& colors, const GridF& alpha, const std::string& overlayName = "default", int layer = 0);
    AbstractPlotter& setOverlay(const std::pair<GridF, GridF>& colorAndAlpha, const std::string& overlayName = "default", int layer = 0);
    AbstractPlotter& showOverlay(const std::string& overlayName = "default");
    AbstractPlotter& hideOverlay(const std::string& overlayName = "default");

    int exec();
    AbstractPlotter& saveFig(const std::string& filename);
    AbstractPlotter& copyToClipboard();
    void resizeEvent(QResizeEvent* event);
    void showEvent(QShowEvent* event);
    void hideEvent(QHideEvent *event);

    QTimer *animate(std::function<bool()> callback, int interval_ms = 30);

    // AbstractPlotter& reset();

    bool hasPlotValues() const { return !this->dataModel->plotLineData.empty(); }
    bool hasScatterValues() const { return !this->dataModel->scatterData.empty(); }
    bool hasImage() const { return !this->dataModel->getImage().empty(); }
    bool hasVectorField() const { return !this->dataModel->vectorData.field.empty(); }

    // void setOnMousePressed(const std::function<void(const Vector3&, Vector3, bool, bool)>& callback) { this->onMousePressedCallbacks.push_back(callback); }
    // void setOnMouseMoved(const std::function<void(const Vector3&, const Vector3&, QMouseEvent*)>& callback) { this->onMouseMovedCallbacks.push_back(callback); }

    ChartView* chartView;
    std::shared_ptr<PlotModel> dataModel;
    QLabel* mouseInfoLabel = nullptr;

    InterfaceUI* mainInterface;
    InterfaceUI* toolsInterface;
    InterfaceUI* viewOptionsInterface;
    InterfaceUI* saveCopyInterface;
    InterfaceUI* infosInterface;
    InterfaceUI* viewAndCopyInterface = nullptr;

    std::string name;

    virtual AbstractPlotter& updateUI(bool forceUpdate = false);
    virtual void displayInfoUnderMouse(const Vector3& relativeMousePos);

    virtual AbstractPlotter& updateToolsInterface();
    virtual AbstractPlotter& updateViewOptionsInterface();
    virtual AbstractPlotter& updateSaveCopyInterface();
    virtual AbstractPlotter& updateInfosInterface();

protected:
    static std::string defaultName;
    static std::map<std::string, AbstractPlotter*> instances;
    static std::string id_name;

    DECLARE_EVENT(MouseMoved, (const Vector3& mousePos, const Vector3& previousPos, QMouseEvent* event), (mousePos, previousPos, event))
    DECLARE_EVENT(MousePressed, (const Vector3& mousePos, Vector3 value, bool leftClick, bool rightClick), (mousePos, value, leftClick, rightClick))
    DECLARE_EVENT(MouseReleased, (const Vector3& mousePos), (mousePos))
    DECLARE_EVENT(Update, (), ())
    DECLARE_EVENT(UserModifiedData, (const std::vector<size_t>& affectedSeries), (affectedSeries))

public Q_SLOTS:
    virtual void draw();
    virtual void show();
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
