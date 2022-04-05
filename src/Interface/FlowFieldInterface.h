#ifndef FLOWFIELDINTERFACE_H
#define FLOWFIELDINTERFACE_H

class FlowFieldInterface;
#include <QWidget>
#include "Interface/CustomInteractiveObject.h"
#include "TerrainGen/VoxelGrid.h"

class FlowFieldInterface : public CustomInteractiveObject
{
    Q_OBJECT
public:
    FlowFieldInterface();

    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);

    QLayout* createGUI();

public Q_SLOTS:
    void show();
    void hide();

public:
    std::shared_ptr<VoxelGrid> voxelGrid;

protected:
    QHBoxLayout* flowFieldLayout;
    QPushButton* flowFieldComputeButton;
    QCheckBox* flowFieldDisplayButton;
};

#endif // FLOWFIELDINTERFACE_H
