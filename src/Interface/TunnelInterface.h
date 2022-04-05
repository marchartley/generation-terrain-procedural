#ifndef TUNNELINTERFACE_H
#define TUNNELINTERFACE_H

class TunnelInterface;
#include "Interface/CustomInteractiveObject.h"
#include "TerrainGen/VoxelGrid.h"
#include "Interface/FancySlider.h"

class TunnelInterface : public CustomInteractiveObject
{
    Q_OBJECT
public:
    TunnelInterface();

    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);

    QLayout* createGUI();

public Q_SLOTS:
    void show();
    void hide();

public:
    std::shared_ptr<VoxelGrid> voxelGrid;

protected:
    QHBoxLayout* tunnelLayout;
    QPushButton* addControlPointButton;
    QPushButton* tunnelClearControlPointButton;
    FancySlider* tunnelSizeSlider;
    FancySlider* tunnelStrengthSlider;
    QPushButton* tunnelCreateMatter;
    QPushButton* tunnelRemoveMatter;
    QPushButton* tunnelCreateCrack;
    QCheckBox* tunnelDisplayButton;
};

#endif // TUNNELINTERFACE_H
