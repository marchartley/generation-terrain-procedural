#ifndef UNDOREDOINTERFACE_H
#define UNDOREDOINTERFACE_H

#include "Interface/ActionInterface.h"
#include "TerrainGen/VoxelGrid.h"

class UndoRedoInterface : public ActionInterface
{
    Q_OBJECT
public:
    UndoRedoInterface(QWidget *parent = nullptr);

    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);

    void display();
    void replay(nlohmann::json action);

    void keyPressEvent(QKeyEvent* event);

    QLayout* createGUI();

Q_SIGNALS:
    void updated();

public Q_SLOTS:
    void show();
    void hide();

    bool undo();
    bool redo();

public:
    std::shared_ptr<VoxelGrid> voxelGrid;

protected:
};

#endif // UNDOREDOINTERFACE_H
