#ifndef FAULTSLIPINTERFACE_H
#define FAULTSLIPINTERFACE_H

class FaultSlipInterface;
#include "Interface/ControlPoint.h"
#include "Interface/InteractiveVector.h"
#include "Utils/BSpline.h"
#include <QWidget>
#include "TerrainGen/VoxelGrid.h"
#include "TerrainModification/FaultSlip.h"
#include "Interface/ActionInterface.h"

class FaultSlipInterface : public ActionInterface
{
    Q_OBJECT
public:
    FaultSlipInterface(QWidget *parent = nullptr);

    void display();
    void remesh();

    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);

    void replay(nlohmann::json action);
    std::shared_ptr<VoxelGrid> voxelGrid;
//    std::shared_ptr<FaultSlip> faultSlip;
    FaultSlip faultSlip;

//    std::vector<BSpline> karstPaths;

    Mesh planeMesh;

    QLayout* createGUI();

Q_SIGNALS:
    void faultSlipApplied();

public Q_SLOTS:
    void updateSlipVector(Vector3 newSlipVector = Vector3());
    void updatePoints();
    void computeFaultSlip();
    void setSideAffected(bool isRightSide);

    void hide();
    void show();

protected:
    std::unique_ptr<ControlPoint> firstSlipControlPoint;
//    ControlPoint *secondSlipControlPoint;
    std::unique_ptr<InteractiveVector> slipVector;

    QHBoxLayout* faultSlipLayout = nullptr;
    QPushButton* faultApplyButton;
    QCheckBox* faultSideApplied;
    QCheckBox* faultDisplayButton;

    void setBindings();
};

#endif // FAULTSLIPINTERFACE_H
