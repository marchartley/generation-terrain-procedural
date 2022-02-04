#ifndef INTERFACE_H
#define INTERFACE_H

#include <QObject>
#include <QSlider>
#include <QGridLayout>
#include "Interface/Viewer.h"
#include "Interface/RangeSlider.h"
#include "Interface/FancySlider.h"
#include "Interface/Spoiler.h"

#include "Interface/KarstPathGenerationInterface.h"
#include "Interface/SpaceColonizationInterface.h"

class ViewerInterface : public QDialog{
public:
    ViewerInterface();
    ~ViewerInterface();

    void setupUi(QDialog *Dialog);

    void setupBindings(QDialog* Dialog);

    void retranslateUi(QDialog *Dialog);

    void setAllValuesToFitViewerDefaults(Viewer* viewer);

    void closeEvent(QCloseEvent* e);

    QGridLayout* mainLayout;
    QVBoxLayout* controlLayout;
    Viewer* viewer;

    std::shared_ptr<KarstPathGenerationInterface> karstPathGeneration;
    std::shared_ptr<SpaceColonizationInterface> spaceColonization;


    QHBoxLayout* loadSaveLayout;
    Spoiler* loadSaveBox;
    QHBoxLayout* randomRocksLayout;
    Spoiler* randomRocksBox;
    QHBoxLayout* currentSimulationLayout;
    Spoiler* currentSimulationBox;
    QHBoxLayout* manualRocksLayout;
    Spoiler* manualRocksBox;
    QHBoxLayout* curvesErosionLayout;
    Spoiler* curvesErosionBox;
    QHBoxLayout* karstCreationLayout;
    Spoiler* karstCreationBox;
    QHBoxLayout* spaceColonizerLayout;
    Spoiler* spaceColonizerBox;
    QHBoxLayout* gravityLayout;
    Spoiler* gravityBox;
    QHBoxLayout* recordingLayout;
    Spoiler* recordingBox;
    QHBoxLayout* displayModeLayout;
    Spoiler* displayModeBox;
    QHBoxLayout* algorithmLayout;
    Spoiler* algorithmBox;
    QHBoxLayout* displayTypeLayout;
    Spoiler* displayTypeBox;
    QHBoxLayout* LoDChooserLayout;
    Spoiler* LoDChooserBox;
    QPushButton* loadButton;
    QPushButton* saveButton;
    FancySlider* randomRocksSizeSlider;
    FancySlider* randomRocksStrengthSlider;
    FancySlider* randomRocksQuantitySlider;
    QPushButton* sendRandomRocksFromCam;
    QPushButton* sendRandomRocks;
    QCheckBox* displayRandomRocks;
    QCheckBox* displayFailedRandomRocks;
    FancySlider* currentSimulationFlowfieldStrengthSlider;
    FancySlider* currentSimulationStrengthSlider;
    FancySlider* currentSimulationRandomDirectionSlider;
    QCheckBox* displayFlowfield;
    QPushButton* recomputeFlowfieldButton;
    FancySlider* manualRocksSizeSlider;
    FancySlider* manualRocksStrengthSlider;
    QRadioButton* addingMode;
    QRadioButton* suppressMode;
    QPushButton* curvesAddControlPointButton;
    QPushButton* curvesClearControlPointButton;
    FancySlider* curvesErosionSizeSlider;
    FancySlider* curvesErosionStrengthSlider;
    QPushButton* curvesErosionCreateMatter;
    QPushButton* curvesErosionRemoveMatter;
    QCheckBox* displayCurvesErosion;

    QPushButton* gravityGlobalButton;
    QPushButton* gravitySandButton;
    QPushButton* startStopRecording;
    QRadioButton* wireModeButton;
    QRadioButton* fillModeButton;
    QRadioButton* invisibleModeButton;
    RangeSlider* mapSliceSliderX;
    RangeSlider* mapSliceSliderY;
    RangeSlider* mapSliceSliderZ;
    QRadioButton* noAlgorithmButton;
    QRadioButton* marchingCubesButton;
    QRadioButton* gridModeButton;
    QRadioButton* voxelsModeButton;
    QRadioButton* layerModeButton;
    FancySlider* LoDChooserSlider;
    QPushButton* LoDChooserConfirmButton;

    std::map<QWidget*, QWidget*> widgetsMapping;
};

#endif // INTERFACE_H
