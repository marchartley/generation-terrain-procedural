#ifndef ENVOBJSINTERFACE_H
#define ENVOBJSINTERFACE_H

#include <QWidget>
#include "Interface/ActionInterface.h"

class EnvObjsInterface : public ActionInterface
{
    Q_OBJECT
public:
    EnvObjsInterface(QWidget *parent = nullptr);

    void affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch = nullptr);

    void display(const Vector3& camPos = Vector3(false));

    void replay(nlohmann::json action);

    QLayout* createGUI();

public Q_SLOTS:
    void show();
    void hide();
    virtual void afterTerrainUpdated();

    void instantiateObject();

public:

};


#endif // ENVOBJSINTERFACE_H
