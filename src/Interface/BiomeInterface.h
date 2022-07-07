#ifndef BIOMEINTERFACE_H
#define BIOMEINTERFACE_H


#include "Biomes/BiomeInstance.h"
class BiomeInterface;
#include "Interface/ControlPoint.h"
#include "Interface/InteractiveVector.h"
#include "Utils/BSpline.h"
#include <QWidget>
#include "TerrainGen/VoxelGrid.h"
#include "Interface/ActionInterface.h"

class BiomeReplacementDialog;
class BiomeInterface : public ActionInterface
{
    Q_OBJECT
public:
    BiomeInterface(QWidget *parent = nullptr);

    void display();

    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);
    void affectHeightmap(std::shared_ptr<Grid> heightmap);

    void replay(nlohmann::json action);

    void mouseMoveEvent(QMouseEvent* event);
    void keyPressEvent(QKeyEvent* event);
    void keyReleaseEvent(QKeyEvent* event);
    void wheelEvent(QWheelEvent* event);
    void mousePressEvent(QMouseEvent* event);

    std::shared_ptr<VoxelGrid> voxelGrid;
    std::shared_ptr<Grid> heightmap;

    QLayout* createGUI();

    std::vector<Mesh> selectionPlanes;
    BiomeModel biomeModel;
    int tempIndex = -1;

Q_SIGNALS:

public Q_SLOTS:

    void generateBiomes(std::shared_ptr<BiomeInstance> predefinedBiomeInstance = nullptr);
    void randomize();

    void replaceBiome(std::shared_ptr<BiomeInstance> biomeToReplace, BiomeInstance newBiome);

    void hide();
    void show();

    void mouseClickedOnMapEvent(Vector3 mousePosInMap, bool mouseInMap, QMouseEvent* event);
    void updateSelectionPlaneToFitBiome(int biomeID, int planeIndex);

    void displayUniqueSelection(int selectionIndex);

protected:

    void setBindings();
    void updateBiomeSelectionGui();

    std::vector<int> selectedBiomeIDs;
//    int selectedBiomeID = -1;
    QLayout* layout = nullptr;

    QListWidget* biomeSelectionGui = nullptr;
//    QLayout* biomeSelectionGuiLayout = nullptr;
public:
    std::shared_ptr<BiomeInstance> rootBiome;
    std::vector<BiomeInstance> possibleBiomeInstances;

    BiomeReplacementDialog* replaceDialog;
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
    QListWidget* allAvailableBiomes;
    QPushButton* cancelButton;
    QPushButton* validButton;
    BiomeInterface* caller = nullptr;
    int selectedBiomeIndex = -1;
};

#endif // BIOMEINTERFACE_H
