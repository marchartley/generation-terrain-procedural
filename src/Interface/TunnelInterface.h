#ifndef TUNNELINTERFACE_H
#define TUNNELINTERFACE_H

class TunnelInterface;
#include "Interface/ActionInterface.h"
#include "TerrainGen/VoxelGrid.h"
#include "Interface/FancySlider.h"
#include "Interface/ControlPoint.h"
#include "Karst/KarstHole.h"

class TunnelInterface : public ActionInterface
{
    Q_OBJECT
public:
    TunnelInterface(QWidget *parent = nullptr);

    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);

    void display();
    void replay(nlohmann::json action);


    QLayout* createGUI();

Q_SIGNALS:
    void needToClipView(Vector3 direction, Vector3 center, bool active);

public Q_SLOTS:
    void show();
    void hide();

//    void setTunnelSize(int newSize) { this->tunnelSize = newSize; computeTunnelPreview(); }
    void setTunnelWidth(int newSize);
    void setTunnelHeight(int newSize);
    void setErosionStrength(float newStrength) { this->erosionStrength = newStrength;}
    void addCurvesControlPoint(Vector3 pos, bool justUpdatePath = false);

    void updateStartingShape();
    void updateEndingShape();

    void clearTunnelPoints();
    void createTunnel(bool removingMatter = true);
    void createCrack(bool removingMatter = true);

    void mouseClickInWorldEvent(Vector3 mousePosInWorld, bool mouseInMap, QMouseEvent *event);
    void wheelEvent(QWheelEvent* event);

public:
    std::shared_ptr<VoxelGrid> voxelGrid;

    KarstHolePredefinedShapes startingShape;
    KarstHolePredefinedShapes endingShape;

protected:
    float tunnelWidth = 10, tunnelHeight = 10;
    float erosionStrength = 3.f;

    void computeTunnelPreview();

    std::vector<Vector3> currentTunnelPoints;
    std::vector<std::unique_ptr<ControlPoint>> controlPoints;
    Mesh tunnelPreview;

    QHBoxLayout* tunnelLayout;
    QPushButton* addControlPointButton;
    QPushButton* tunnelClearControlPointButton;
    FancySlider* tunnelWidthSlider;
    FancySlider* tunnelHeightSlider;
    FancySlider* tunnelStrengthSlider;
    QPushButton* tunnelCreateMatter;
    QPushButton* tunnelRemoveMatter;
    QPushButton* tunnelCreateCrack;
    QCheckBox* tunnelDisplayButton;

    QComboBox* startingShapeCombobox;
    QComboBox* endingShapeCombobox;
};

#endif // TUNNELINTERFACE_H
