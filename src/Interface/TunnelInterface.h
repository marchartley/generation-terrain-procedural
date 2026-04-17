#ifndef TUNNELINTERFACE_H
#define TUNNELINTERFACE_H

class TunnelInterface;
#include "Interface/ActionInterface.h"
// #include "TerrainGen/VoxelGrid.h"
// #include "GUIElements/FancySlider.h"
#include "GUI3DElements/ControlPoint.h"
#include "Karst/KarstHole.h"

class TunnelInterface : public ActionInterface
{
    Q_OBJECT
public:
    TunnelInterface(QWidget *parent = nullptr);

//    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);

    void display(const Vector3& camPos = Vector3::invalid);
    void replay(nlohmann::json action);


    InterfaceUI* createGUI();
    virtual void affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch = nullptr);

Q_SIGNALS:
    void needToClipView(const Vector3& direction, const Vector3& center, bool active);
    void tunnelCreated(KarstHole& tunnel);

public Q_SLOTS:
    void show();
    void hide();

//    void setTunnelSize(int newSize) { this->tunnelSize = newSize; computeTunnelPreview(); }
    void setTunnelWidth(int newSize);
    void setTunnelHeight(int newSize);
    void setErosionStrength(float newStrength) { this->erosionStrength = newStrength;}
    void addCurvesControlPoint(const Vector3& pos, bool justUpdatePath = false);

    void updateStartingShape();
    void updateEndingShape();

    void clearTunnelPoints();
    void createTunnel(bool removingMatter = true);
    void createCrack(bool removingMatter = true);

    void mouseClickedOnMapEvent(const Vector3& mousePosInWorld, bool mouseInMap, QMouseEvent *event, TerrainModel *model);
    void wheelEvent(QWheelEvent* event);

public:

    KarstHolePredefinedShapes startingShape;
    KarstHolePredefinedShapes endingShape;

protected:
    float tunnelWidth = 10, tunnelHeight = 10;
    float erosionStrength = 3.f;

    void computeTunnelPreview();

    std::vector<Vector3> currentTunnelPoints;
    std::vector<std::shared_ptr<ControlPoint>> controlPoints;
    Mesh tunnelPreview;

    std::vector<ComboboxLineElement<KarstHolePredefinedShapes>*> shapes;
    int startingShapeIndex;
    int endingShapeIndex;
};

#endif // TUNNELINTERFACE_H
