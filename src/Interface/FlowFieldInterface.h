#ifndef FLOWFIELDINTERFACE_H
#define FLOWFIELDINTERFACE_H

class FlowFieldInterface;
#include <QWidget>
#include "Interface/ActionInterface.h"
#include "TerrainGen/VoxelGrid.h"

class FlowFieldInterface : public ActionInterface
{
    Q_OBJECT
public:
    FlowFieldInterface(QWidget *parent = nullptr);

    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);

    void display();

    void replay(nlohmann::json action);

    QLayout* createGUI();

Q_SIGNALS:
    void updated();

public Q_SLOTS:
    void show();
    void hide();

    void recomputeFlowfield();
    void updateFlowfieldDebugMesh();

public:
    std::shared_ptr<VoxelGrid> voxelGrid;

protected:
    Mesh flowMesh;

    QHBoxLayout* flowFieldLayout;
    QPushButton* flowFieldComputeButton;
    QCheckBox* flowFieldDisplayButton;
};

#endif // FLOWFIELDINTERFACE_H
