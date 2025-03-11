#ifndef INTERFACE_H
#define INTERFACE_H

#include <QObject>
#include <QSlider>
#include <QGridLayout>
#include "Interface/Viewer.h"
#include "Interface/RangeSlider.h"
#include "Interface/FancySlider.h"
#include "Interface/Spoiler.h"

#include "Interface/StickyFrame.h"

class ViewerInterface : public QMainWindow{
public:
    ViewerInterface(std::string preloaded_heightmap = "", MapMode displayMode = MapMode::GRID_MODE);
    virtual ~ViewerInterface();

    void setupUi();

    void setupBindings();

    void retranslateUi();

    void setAllValuesToFitViewerDefaults(Viewer* viewer);

    void closeEvent(QCloseEvent* e);

    void resizeEvent(QResizeEvent* e);

    void focusOutEvent(QFocusEvent* e);

public Q_SLOTS:
//    void openInterface(std::string interfaceName, std::shared_ptr<ActionInterface> object);
    void openInterface(std::shared_ptr<ActionInterface> object);

    void hideAllInteractiveParts();

    void erosionsTests();

//    void openMapUI();
//    void saveMapUI();

public:
    StickyFrame* frame;
    std::string lastPanelOpenedByStickyFrame;

    QGridLayout* mainLayout;
    Viewer* viewer;

    /*
    std::shared_ptr<KarstPathGenerationInterface> karstPathGeneration;
    std::shared_ptr<SpaceColonizationInterface> spaceColonization;
    std::shared_ptr<FaultSlipInterface> faultSlip;
    std::shared_ptr<FlowFieldInterface> flowField;
    std::shared_ptr<TunnelInterface> tunnelInterface;
    std::shared_ptr<ManualEditionInterface> manualEditionInterface;
    std::shared_ptr<GravityInterface> gravityInterface;
    std::shared_ptr<UndoRedoInterface> undoRedoInterface;
    std::shared_ptr<TerrainGenerationInterface> terrainGenerationInterface;
    std::shared_ptr<ErosionInterface> erosionInterface;
    std::shared_ptr<HeightmapErosionInterface> heightmapErosionInterface;
    std::shared_ptr<BiomeInterface> biomeInterface;
    std::shared_ptr<SmoothInterface> smoothInterface;
    std::shared_ptr<PrimitivePatchesInterface> patchesInterface;
    std::shared_ptr<TerrainSavingInterface> savingInterface;
    std::shared_ptr<MeshInstanceAmplificationInterface> meshInstanceAmplificationInterface;
    std::shared_ptr<SPHSimulationInterface> sphSimulationInterface;
    std::shared_ptr<FLIPSimulationInterface> flipSimulationInterface;
    std::shared_ptr<WarpFluidSimulationInterface> warpFluidSimulationInterface;
    std::shared_ptr<LBMFluidSimulationInterface> LbmFluidSimulationInterface;
    std::shared_ptr<CoralIslandGeneratorInterface> coralIslandGeneratorInterface;
    std::shared_ptr<SpheroidalErosionInterface> spheroidalErosionInterface;
    */

    std::map<std::string, std::shared_ptr<ActionInterface>> actionInterfaces;

    QMenuBar* toolbox;

    QHBoxLayout* displayModeLayout;
    Spoiler* displayModeBox;
    QHBoxLayout* LoDChooserLayout;
    Spoiler* LoDChooserBox;
    RangeSlider* mapSliceSliderX;
    RangeSlider* mapSliceSliderY;
    RangeSlider* mapSliceSliderZ;
    RangeSlider* isolevelSelectionSlider;
    QCheckBox* mapSliceSmooth;

    std::shared_ptr<ActionInterface> lastOpenedInterface;

//    std::string mapSavingFolder = "saved_maps/";
    std::shared_ptr<std::vector<nlohmann::json>> actionsOnMap;

    float randomValue = 10.f;
};

#endif // INTERFACE_H
