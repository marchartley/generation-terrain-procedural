#ifndef TUNNELINTERFACE_H
#define TUNNELINTERFACE_H

class TunnelInterface;
#include "Interface/CustomInteractiveObject.h"
#include "TerrainGen/VoxelGrid.h"
#include "Interface/FancySlider.h"
#include "Interface/ControlPoint.h"

class TunnelInterface : public CustomInteractiveObject
{
    Q_OBJECT
public:
    TunnelInterface(QWidget *parent = nullptr);

    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);

    void display();


    QLayout* createGUI();

Q_SIGNALS:
    void updated();

public Q_SLOTS:
    void show();
    void hide();

    void setTunnelSize(int newSize) { this->tunnelSize = newSize; computeTunnelPreview(); }
    void setErosionStrength(float newStrength) { this->erosionStrength = newStrength;}
    void addCurvesControlPoint(Vector3 pos, bool justUpdatePath = false);

    void clearTunnelPoints();
    void createTunnel(bool removingMatter = true);
    void createCrack(bool removingMatter = true);

    void mouseClickInWorldEvent(Vector3 mousePosInWorld, bool mouseInMap, QMouseEvent *event);

public:
    std::shared_ptr<VoxelGrid> voxelGrid;

protected:
    float tunnelSize = 10;
    float erosionStrength = 3.f;

    void computeTunnelPreview();

    std::vector<Vector3> currentTunnelPoints;
    std::vector<std::unique_ptr<ControlPoint>> controlPoints;
    Mesh tunnelPreview;

    QHBoxLayout* tunnelLayout;
    QPushButton* addControlPointButton;
    QPushButton* tunnelClearControlPointButton;
    FancySlider* tunnelSizeSlider;
    FancySlider* tunnelStrengthSlider;
    QPushButton* tunnelCreateMatter;
    QPushButton* tunnelRemoveMatter;
    QPushButton* tunnelCreateCrack;
    QCheckBox* tunnelDisplayButton;
};

#endif // TUNNELINTERFACE_H
