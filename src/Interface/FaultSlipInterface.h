#ifndef FAULTSLIPINTERFACE_H
#define FAULTSLIPINTERFACE_H

class FaultSlipInterface;
#include "GUI3DElements/ControlPoint.h"
#include "GUI3DElements/InteractiveVector.h"
#include "Utils/BSpline.h"
// #include <QWidget>
#include "TerrainGen/VoxelGrid.h"
#include "TerrainModification/FaultSlip.h"
#include "Interface/ActionInterface.h"

class FaultSlipInterface : public ActionInterface
{
    Q_OBJECT
public:
    FaultSlipInterface(QWidget *parent = nullptr);

    void display(const Vector3& camPos = Vector3::invalid);
    void remesh();

    void affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch = nullptr);

    void replay(nlohmann::json action);
//    std::shared_ptr<VoxelGrid> voxelGrid;
//    std::shared_ptr<FaultSlip> faultSlip;
    FaultSlip faultSlip;

//    std::vector<BSpline> karstPaths;

    Mesh planeMesh;

    InterfaceUI* createGUI();

Q_SIGNALS:
//    void faultSlipApplied();

public Q_SLOTS:
    void updateSlipVector(const Vector3& newSlipVector = Vector3());
    void updatePoints();
    void computeFaultSlip();
    void setSideAffected(bool isRightSide);

    void hide();
    void show();

protected:
    std::shared_ptr<ControlPoint> firstSlipControlPoint;
//    ControlPoint3D *secondSlipControlPoint;
    std::shared_ptr<InteractiveVector> slipVector;

//    QHBoxLayout* faultSlipLayout = nullptr;
//    QPushButton* faultApplyButton;
//    QCheckBox* faultSideApplied;
//    QCheckBox* faultDisplayButton;

    void setBindings();
};

#endif // FAULTSLIPINTERFACE_H
