#ifndef PLOTTER_H
#define PLOTTER_H

#include "Graphics/AbstractPlotter.h"

class Plotter : public AbstractPlotter {
    Q_OBJECT
private: // Singleton
    Plotter(std::string name, QWidget* parent = nullptr);
    Plotter(std::string name, ChartView* chartView, QWidget* parent = nullptr);

public:
    static Plotter* getInstance(std::string name = "");
    static Plotter* get(std::string name = "") { return Plotter::getInstance(toLower(name)); }
    static Plotter* init(std::string name, ChartView* chartView = nullptr, QWidget* parent = nullptr);

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
