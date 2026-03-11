#ifndef UNDOREDOINTERFACE_H
#define UNDOREDOINTERFACE_H

#include "Interface/ActionInterface.h"
// #include "TerrainGen/VoxelGrid.h"

class UndoRedoInterface : public ActionInterface
{
    Q_OBJECT
public:
    UndoRedoInterface(QWidget *parent = nullptr);

    void display(const Vector3& camPos = Vector3::invalid);
    void replay(nlohmann::json action);

    void keyPressEvent(QKeyEvent* event);

    InterfaceUI* createGUI();

public Q_SLOTS:
    void show();
    void hide();

    bool undo();
    bool redo();

protected:
};

#endif // UNDOREDOINTERFACE_H
