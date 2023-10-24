#ifndef SMOOTHINTERFACE_H
#define SMOOTHINTERFACE_H

#include "ActionInterface.h"

class SmoothInterface : public ActionInterface
{
    Q_OBJECT
public:
    SmoothInterface(QWidget *parent = nullptr);

//    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);

    void display(const Vector3& camPos = Vector3(false));
    void replay(nlohmann::json action);

    bool applySmooth();

    QLayout* createGUI();

public Q_SLOTS:
    void show();
    void hide();

public:
//    std::shared_ptr<VoxelGrid> voxelGrid;

protected:
    QHBoxLayout* smoothLayout;
    QPushButton* smoothComputeButton;
};

#endif // SMOOTHINTERFACE_H
