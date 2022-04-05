#ifndef MANUALEDITIONINTERFACE_H
#define MANUALEDITIONINTERFACE_H

class ManualEditionInterface;

#include "Interface/CustomInteractiveObject.h"
#include <QObject>
#include "TerrainGen/VoxelGrid.h"
#include "Interface/FancySlider.h"

class ManualEditionInterface : public CustomInteractiveObject
{
    Q_OBJECT
public:
    ManualEditionInterface();

    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);
    QLayout* createGUI();

public Q_SLOTS:
    void show();
    void hide();

public:
    std::shared_ptr<VoxelGrid> voxelGrid;

    QHBoxLayout* manualEditLayout;
    FancySlider* manualEditSizeSlider;
    FancySlider* manualEditStrengthSlider;
    QRadioButton* addingModeButton;
    QRadioButton* suppressModeButton;
};

#endif // MANUALEDITIONINTERFACE_H
