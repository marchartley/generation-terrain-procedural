#ifndef GANPAINTERINTERFACE_H
#define GANPAINTERINTERFACE_H

#include "Interface/Interface.h"

class GANPainterInterface : public ActionInterface
{
    Q_OBJECT
public:
    GANPainterInterface(QWidget *parent = nullptr);
    ~GANPainterInterface();

    void display(const Vector3& camPos = Vector3::invalid);
    void replay(nlohmann::json action);

    void affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch = nullptr);

    void runGANs();

    GridV3 computeNewHeightmap();

    GridV3 randomIslandGenerator();

    InterfaceUI* createGUI();

    std::map<std::string, QProcess*> ganProcesses;
    std::map<std::string, std::string> ganModelNames;
    std::map<std::string, float> ganModelWeights;
    std::map<std::string, GridF> ganComputedHeightfields;


    float subsidence = .5f;
};

#endif // GANPAINTERINTERFACE_H
