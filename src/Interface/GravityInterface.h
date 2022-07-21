#ifndef GRAVITYINTERFACE_H
#define GRAVITYINTERFACE_H

class GravityInterface;
#include <QWidget>
#include "Interface/ActionInterface.h"
#include "TerrainGen/VoxelGrid.h"

class GravityInterface : public ActionInterface
{
    Q_OBJECT
public:
    GravityInterface(QWidget *parent = nullptr);

    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);

    void display();
    void replay(nlohmann::json action);

    bool createGlobalGravity();
    bool createSandGravity();

    QLayout* createGUI();

public Q_SLOTS:
    void show();
    void hide();

public:
    std::shared_ptr<VoxelGrid> voxelGrid;

protected:
    QHBoxLayout* gravityLayout;
    QPushButton* gravityComputeButton;
//    QCheckBox* gravityDisplayButton;
};

#endif // GRAVITYINTERFACE_H
