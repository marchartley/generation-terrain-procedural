#include "UndoRedoInterface.h"

#include <QKeyEvent>

UndoRedoInterface::UndoRedoInterface(QWidget *parent) : ActionInterface("undo-redo", parent)
{
    this->createGUI();
}

void UndoRedoInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;
}

void UndoRedoInterface::display()
{
    // Nothing to display
}

void UndoRedoInterface::replay(nlohmann::json action)
{
    if (this->isConcerned(action)) {
        auto& parameters = action.at("parameters");
        bool applyUndo = parameters.at("undo").get<bool>();

        if (applyUndo) {
            this->voxelGrid->undo();
        } else {
            this->voxelGrid->redo();
        }
    }
}

void UndoRedoInterface::keyPressEvent(QKeyEvent *event)
{
    if (event->modifiers().testFlag(Qt::ControlModifier)) {
        if (event->key() == Qt::Key_Z)
            undo();
        if (event->key() == Qt::Key_Y)
            redo();
    }
}


bool UndoRedoInterface::undo()
{
    this->voxelGrid->undo();

    this->addTerrainAction(nlohmann::json({
                                               {"undo", true},
                                            }));
    Q_EMIT updated();
    return false;
}

bool UndoRedoInterface::redo()
{
    this->voxelGrid->redo();
    Q_EMIT updated();


    this->addTerrainAction(nlohmann::json({
                                           {"undo", false},
                                          }));
    return false;
}

void UndoRedoInterface::hide()
{
    CustomInteractiveObject::hide();
}

void UndoRedoInterface::show()
{
    CustomInteractiveObject::show();
}

QLayout* UndoRedoInterface::createGUI()
{
    QHBoxLayout* nothing = new QHBoxLayout;
    return nothing;
}
