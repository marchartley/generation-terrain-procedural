#include "ManualEditionInterface.h"

#include "GUIElements/InterfaceUtils.h"
#include "TerrainModification/RockErosion.h"

// #include "serialization/Serializer.h"

ManualEditionInterface::ManualEditionInterface(QWidget *parent)
    : ActionInterface("manualedit", "Manual edit", "digging", "Manual editing", "manual-edit_button.png", parent)
{
    this->grabber = std::make_shared<ControlPoint>(Vector3(), this->manualEditionSize/2.f, NEUTRAL, false);
    this->grabber->setGrabberStateColor(POSITIVE, std::vector<float>({.1f, .9f, .1f, .2f}));
    this->grabber->setGrabberStateColor(NEGATIVE, std::vector<float>({.9f, .1f, .1f, .2f}));
    setAddingMode(addingMode);
}

void ManualEditionInterface::display(const Vector3& camPos)
{
    if (!this->isVisible())
        return;

//    std::cout << (readyToModify ? "modif" : "not modif") << std::endl;
    if (this->readyToModify) {
        this->grabber->setState((this->addingMode ? POSITIVE : NEGATIVE));
    } else {
        this->grabber->setState(NEUTRAL);
    }
    grabber->display();
}

void ManualEditionInterface::replay(nlohmann::json action)
{
    if (this->isConcerned(action)) {
        auto& parameters = action.at("parameters");

        float size = parameters.at("size").get<float>() * random_gen::generate(0.f, 2.f);
        float strength = parameters.at("strength").get<float>() * random_gen::generate(0.f, 2.f);
        Vector3 position = parameters.at("position").get<Vector3>() + Vector3::random(0.f, 20.f);
        bool addingMode = parameters.at("addingMode").get<bool>();

        RockErosion rock(size, strength);
        rock.Apply(voxelGrid, position, addingMode);
    }
}

/*void ManualEditionInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;
}*/

void ManualEditionInterface::show()
{
    grabber->show();
    CustomInteractiveObject::show();
}

void ManualEditionInterface::hide()
{
    grabber->hide();
    CustomInteractiveObject::hide();
}

void ManualEditionInterface::setSize(int size)
{
    if (size < 0) size = 0;
    this->manualEditionSize = size;
    this->grabber->setRadius(size/2.f);
}

void ManualEditionInterface::setStrength(float strength)
{
    this->manualEditionStrength = strength;
}
void ManualEditionInterface::setAddingMode(bool newMode)
{
    this->addingMode = newMode;
//    this->grabber->state = (this->addingMode ? POSITIVE : NEGATIVE);
}

void ManualEditionInterface::setPosition(const Vector3& newPosition)
{
    this->grabber->move(newPosition);
}

void ManualEditionInterface::applyModification()
{
    float size = this->manualEditionSize;
    float strength = this->manualEditionStrength;
    Vector3 position = this->grabber->getPosition();

//    RockErosion rock(size, strength);
//    rock.Apply(voxelGrid, position, addingMode);
    GridF modif = RockErosion::createPrecomputedAttackMask(size);
    modif *= strength / (modif.at(size * .5f, size * .5f, size * .5f));
    voxelGrid->applyModification(modif * (addingMode ? 1.f : -1.f), position - Vector3(size * .5f, size * .5f, size * .5f));

    this->addTerrainAction(nlohmann::json({
                                           {"size", size},
                                           {"strength", strength},
                                           {"position", position},
                                           {"addingMode", addingMode}
                                          }));

    Q_EMIT this->updated();
}

void ManualEditionInterface::mouseMovedOnMapEvent(const Vector3& mouseWorldPosition, TerrainModel *model)
{
    this->setPosition(mouseWorldPosition);
    if (readyToModify) {
//        this->applyModification();
    }
}

void ManualEditionInterface::mouseClickedOnMapEvent(const Vector3& mousePosInMap, bool mouseInMap, QMouseEvent *event, TerrainModel* model)
{
    if (this->isVisible() && mouseInMap &&
            event->button() == Qt::MouseButton::LeftButton &&
            readyToModify)
    {
        this->setPosition(mousePosInMap);
        this->applyModification();
    }
}

void ManualEditionInterface::mouseMoveEvent(QMouseEvent *event)
{
    // Just used to cancel the Ctrl and Alt buttons
    if (readyToModify) {
        if (!event->modifiers().testFlag(Qt::ControlModifier) && !event->modifiers().testFlag(Qt::ShiftModifier))
        {
            readyToModify = false;
        } else {
        }
    }
//    CustomInteractiveObject::mouseMoveEvent(event);
}

void ManualEditionInterface::keyPressEvent(QKeyEvent *event)
{
//    std::cout << "Here!" << (event->isAutoRepeat() ? " on auto " : "not auto") << " - " << (event->key() == Qt::Key_Shift || event->key() == Qt::Key_Control ? "good key" : "bad key") << std::endl;
    if (!event->isAutoRepeat()) {
        this->readyToModify = true;
        if (event->key() == Qt::Key_Shift) {
            this->addingMode = false;
        } else if (event->key() == Qt::Key_Control) {
            this->addingMode = true;
        } else {
            this->readyToModify = false;
        }
    }
    CustomInteractiveObject::keyPressEvent(event);
}

void ManualEditionInterface::keyReleaseEvent(QKeyEvent *event)
{

    if (event->key() == Qt::Key_Shift || event->key() == Qt::Key_Control) {
        this->readyToModify = false;
        Q_EMIT this->updated();
    }
    CustomInteractiveObject::keyReleaseEvent(event);
}

void ManualEditionInterface::wheelEvent(QWheelEvent* event)
{
    if (event->modifiers().testFlag(Qt::ControlModifier)) {
        this->setSize(this->manualEditionSize - event->angleDelta().y() / 2);
        Q_EMIT this->updated();
    }
    CustomInteractiveObject::wheelEvent(event);
}


InterfaceUI* ManualEditionInterface::createGUI()
{
    auto UI = new InterfaceUI();
    // Manual rock erosion layout
    auto manualEditSizeSlider = new SliderElement("Taille", 1, 200, 1, this->manualEditionStrength);
    auto manualEditStrengthSlider = new SliderElement("Force", 0.0, 3.0, 0.1, this->manualEditionStrength);
    auto addingModeButton = new RadioButtonElement("Ajouter de la matière", true, addingMode);
    auto suppressModeButton = new RadioButtonElement("Detruire de la matière", false, addingMode);


    UI->add({
        manualEditSizeSlider,
        // manualEditStrengthSlider,
        addingModeButton,
        suppressModeButton
    });

    return UI;
}
