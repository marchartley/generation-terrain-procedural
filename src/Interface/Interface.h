#ifndef INTERFACE_H
#define INTERFACE_H
#include <QMainWindow>

#include "Interface/Viewer.h"


QT_BEGIN_NAMESPACE
class QGridLayout;
class QMenuBar;
QT_END_NAMESPACE

class StickyFrame;

class ViewerInterface : public QMainWindow{
public:
    ViewerInterface(const std::string& preloaded_heightmap = "", MapMode displayMode = MapMode::GRID_MODE);
    virtual ~ViewerInterface();

    void setupUi();
    void setupBindings();
    void retranslateUi();
    void setAllValuesToFitViewerDefaults(Viewer* viewer);
    void closeEvent(QCloseEvent* e);
    void resizeEvent(QResizeEvent* e);
    void focusOutEvent(QFocusEvent* e);

public Q_SLOTS:
    void openInterface(std::shared_ptr<ActionInterface> object);
    void hideAllInteractiveParts();
    void erosionsTests();

public:
    StickyFrame* frame;
    std::string lastPanelOpenedByStickyFrame;

    QGridLayout* mainLayout;
    Viewer* viewer;

    std::map<std::string, std::shared_ptr<ActionInterface>> actionInterfaces;

    QMenuBar* toolbox;

    InterfaceUI* displayModeLayout;
    // Spoiler* displayModeBox;
    // QHBoxLayout* LoDChooserLayout;
    // Spoiler* LoDChooserBox;
    RangeSliderElement* mapSliceSliderX;
    RangeSliderElement* mapSliceSliderY;
    RangeSliderElement* mapSliceSliderZ;
    RangeSliderElement* isolevelSelectionSlider;
    CheckboxElement* mapSliceSmooth;

    std::shared_ptr<ActionInterface> lastOpenedInterface;

    std::shared_ptr<std::vector<nlohmann::json>> actionsOnMap;

    float randomValue = 10.f;
};

#endif // INTERFACE_H
