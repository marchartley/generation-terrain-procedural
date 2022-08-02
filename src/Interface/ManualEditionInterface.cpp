#include "ManualEditionInterface.h"

#include "Interface/InterfaceUtils.h"
#include "TerrainModification/RockErosion.h"

ManualEditionInterface::ManualEditionInterface(QWidget *parent) : ActionInterface("manualEdit", parent)
{
    this->grabber = std::make_unique<ControlPoint>(Vector3(), this->manualEditionSize/2.f, NEUTRAL, false);
    this->grabber->setGrabberStateColor(POSITIVE, std::vector<float>({.1f, .9f, .1f, .2f}));
    this->grabber->setGrabberStateColor(NEGATIVE, std::vector<float>({.9f, .1f, .1f, .2f}));
    setAddingMode(addingMode);
}

void ManualEditionInterface::display()
{
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

        float size = parameters.at("size").get<float>();
        float strength = parameters.at("strength").get<float>();
        Vector3 position = json_to_vec3(parameters.at("position"));
        bool addingMode = parameters.at("addingMode").get<bool>();

        RockErosion rock(size, strength);
        rock.Apply(voxelGrid, position, addingMode, false);
    }
}

void ManualEditionInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;
}

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
    this->grabber->radius = size/2.f;
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

void ManualEditionInterface::setPosition(Vector3 newPosition)
{
    this->grabber->move(newPosition);
}

void ManualEditionInterface::applyModification()
{
    float size = this->manualEditionSize;
    float strength = this->manualEditionStrength;
    Vector3 position = this->grabber->getPosition();

    RockErosion rock(size, strength);
    rock.Apply(voxelGrid, position, addingMode, true);

    this->addTerrainAction(nlohmann::json({
                                           {"size", size},
                                           {"strength", strength},
                                           {"position", vec3_to_json(position)},
                                           {"addingMode", addingMode}
                                          }));

    Q_EMIT this->updated();
}

void ManualEditionInterface::mouseClickedOnMapEvent(Vector3 mousePosInMap, bool mouseInMap, QMouseEvent *event)
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
    CustomInteractiveObject::mouseMoveEvent(event);
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


QLayout *ManualEditionInterface::createGUI()
{
    this->manualEditLayout = new QHBoxLayout();
    // Manual rock erosion layout
    this->manualEditSizeSlider = new FancySlider(Qt::Orientation::Horizontal, 1, 200);
    this->manualEditStrengthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.0, 3.0, 0.1);
    this->addingModeButton = new QRadioButton("Ajouter de la matière");
    this->suppressModeButton = new QRadioButton("Detruire de la matière");
    manualEditLayout->addWidget(createSliderGroup("Taille", manualEditSizeSlider));
//    manualEditLayout->addWidget(createSliderGroup("Force", manualEditStrengthSlider));
    manualEditLayout->addWidget(createVerticalGroup({addingModeButton, suppressModeButton}));

    this->manualEditSizeSlider->setValue(this->manualEditionSize);
    this->manualEditStrengthSlider->setfValue(this->manualEditionStrength);
    this->addingModeButton->setChecked(this->addingMode);
    this->addingModeButton->setChecked(!this->addingMode);

    QObject::connect(manualEditSizeSlider, &FancySlider::valueChanged, this, &ManualEditionInterface::setSize);
    QObject::connect(manualEditStrengthSlider, &FancySlider::floatValueChanged, this, [=](float newStrength) { this->manualEditionStrength = newStrength; });
    QObject::connect(addingModeButton, &QRadioButton::clicked, this, [=](){ this->setAddingMode(true); } );
    QObject::connect(suppressModeButton, &QRadioButton::clicked, this, [=](){ this->setAddingMode(false); } );

    return this->manualEditLayout;
}
