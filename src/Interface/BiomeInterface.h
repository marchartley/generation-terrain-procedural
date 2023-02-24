#ifndef BIOMEINTERFACE_H
#define BIOMEINTERFACE_H


#include "Biomes/BiomeInstance.h"
class BiomeInterface;
//#include "Interface/ControlPoint.h"
//#include "Interface/InteractiveVector.h"
#include "Utils/BSpline.h"
#include <QWidget>
#include "TerrainGen/VoxelGrid.h"
#include "TerrainGen/Grid.h"
#include "Interface/ActionInterface.h"
#include "Karst/KarstHole.h"
#include "Interface/HierarchicalListWidget.h"

class BiomeReplacementDialog;
class BiomeInterface : public ActionInterface
{
    Q_OBJECT
public:
    BiomeInterface(QWidget *parent = nullptr);

    void display();

//    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);
//    void affectHeightmap(std::shared_ptr<Grid> heightmap);

    void replay(nlohmann::json action);

    void mouseMoveEvent(QMouseEvent* event);
    void keyPressEvent(QKeyEvent* event);
    void keyReleaseEvent(QKeyEvent* event);
    void wheelEvent(QWheelEvent* event);
    void mousePressEvent(QMouseEvent* event);

    QLayout* createGUI();

    std::vector<Mesh> selectionPlanes;
    BiomeModel biomeModel;
    BiomeModel modifiedBiomeModel;
    int tempIndex = -1;
    Vector3 tempMousePos;

    void setVoxelGridSizeFactor(float newFactor);

public Q_SLOTS:
    void displayAllBiomes();
    void interchangeBiomes();
    void mouseDoubleClickOnMapEvent(Vector3 mousePosition, bool mouseInMap, QMouseEvent* event, TerrainModel *model);

    void generateBiomes(std::shared_ptr<BiomeInstance> predefinedBiomeInstance = nullptr);
    void randomize();

    void replaceBiome(std::shared_ptr<BiomeInstance> biomeToReplace, std::shared_ptr<BiomeInstance> newBiome);

    void hide();
    void show();

    void addTunnel(KarstHole& hole);

    void mouseClickedOnMapEvent(Vector3 mousePosInMap, bool mouseInMap, QMouseEvent* event, TerrainModel *model);
    void updateSelectionPlaneToFitBiome(int biomeID, int planeIndex);

    void displayUniqueSelection(int selectionIndex);

    Matrix3<float> prepareTrench(std::shared_ptr<BiomeInstance> biome);
    Matrix3<float> prepareCoralWall(std::shared_ptr<BiomeInstance> biome);
    Matrix3<float> prepareArche(std::shared_ptr<BiomeInstance> biome);
    Matrix3<float> preparePatateCorail(std::shared_ptr<BiomeInstance> biome);

Q_SIGNALS:
    void terrainViewModified(Vector3 newOffset, float newScaling);

protected:
    Vector3 fromHeightmapPosToVoxels(Vector3 pos);
    Vector3 fromVoxelsPosToHeightmap(Vector3 pos);
    void setBindings();
    void updateBiomeSelectionGui();
    void deleteSelectedBiomes();

    void deleteBiomeFromID(int ID);

    std::vector<int> selectedBiomeIDs;
//    int selectedBiomeID = -1;
    QLayout* layout = nullptr;

    HierarchicalListWidget* biomeSelectionGui = nullptr;
//    QLayout* biomeSelectionGuiLayout = nullptr;
public:
    std::shared_ptr<BiomeInstance> rootBiome;
    std::vector<BiomeInstance> possibleBiomeInstances;

    BiomeReplacementDialog* replaceDialog;

    float voxelGridScaleFactor = 1.f;
    Vector3 voxelGridOffsetStart;
};


class BiomeReplacementDialog : public QDialog {
    Q_OBJECT
public:
    BiomeReplacementDialog(BiomeInterface *caller = nullptr);


public Q_SLOTS:
    void open();
    void cancel();
    void confirm();

public:
    HierarchicalListWidget* allAvailableBiomes;
    QPushButton* cancelButton;
    QPushButton* validButton;
    BiomeInterface* caller = nullptr;
    int selectedBiomeIndex = -1;
};

#endif // BIOMEINTERFACE_H
