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

//    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);
    void affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch);

    void display();
    void displayPressureDensities();
    void displayFlows();
    void displaySumOfFlows();

    void replay(nlohmann::json action);

    QLayout* createGUI();

public Q_SLOTS:
    void show();
    void hide();

    void recomputeFlowfield();
    void updateFlowfieldDebugMesh();


protected:
//    Mesh flowMesh;

    std::vector<Mesh> flowMeshes;
    Mesh sumOfFlowsMesh;

    Matrix3<float> pressureDensityVoxels;
    Matrix3<Vector3> totalFlow;
    Mesh pressureDensityMesh;

    std::vector<bool> displayedFlowfields;
    bool displayingPressure = false;
    bool displayingSumOfFlows = false;
};

#endif // FLOWFIELDINTERFACE_H
