#ifndef INTERFACE_H
#define INTERFACE_H

#include <QObject>
#include <QSlider>
#include <QGridLayout>
#include "Interface/Viewer.h"
#include "GUIElements/RangeSlider.h"
#include "GUIElements/FancySlider.h"
#include "GUIElements/Spoiler.h"

#include "GUIElements/StickyFrame.h"

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
