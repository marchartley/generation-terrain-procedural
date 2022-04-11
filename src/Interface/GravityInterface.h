#ifndef GRAVITYINTERFACE_H
#define GRAVITYINTERFACE_H

class GravityInterface;
#include <QWidget>
#include "Interface/CustomInteractiveObject.h"
#include "TerrainGen/VoxelGrid.h"

class GravityInterface : public CustomInteractiveObject
{
    Q_OBJECT
public:
    GravityInterface(QWidget *parent = nullptr);

    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);

    void display();

    bool createGlobalGravity();
    bool createSandGravity();

    QLayout* createGUI();

Q_SIGNALS:
    void updated();

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
