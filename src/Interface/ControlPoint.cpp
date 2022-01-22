#include "ControlPoint.h"

std::shared_ptr<Shader> ControlPoint::base_shader = nullptr;
std::map<GrabberState, std::vector<float>> ControlPoint::GrabberStateColor = {
    {GrabberState::HIDDEN, {.0f, .0f, .0f, 1.f}},
    {GrabberState::INACTIVE, {.3f, .0f, .0f, 1.f}},
    {GrabberState::ACTIVE, {.8f, .0f, .0f, 1.f}}
};

ControlPoint::ControlPoint()
    : ControlPoint(Vector3())
{
    this->isDisplayed = false;
}

ControlPoint::ControlPoint(Vector3 pos, float radius, GrabberState state, bool useTheManipulatedFrame)
    : Mesh((ControlPoint::base_shader ? std::make_shared<Shader>(*ControlPoint::base_shader) : nullptr), true),
      position(pos), radius(radius), state(state), useManipFrame(useTheManipulatedFrame), shape(radius, position, 10, 10)
{
    this->move(pos);
    if (!useManipFrame) {
        this->manipFrame.removeFromMouseGrabberPool();
    } else {
        this->manipFrame.addInMouseGrabberPool();
    }
//    QObject::connect(&this->manipFrame, SIGNAL(manipulated()), this, SLOT(updateSphere()));
//    std::cout << this->manipFrame.isInMouseGrabberPool() << " " << this->manipFrame.grabsMouse() << std::endl;
}

ControlPoint::~ControlPoint()
{
    this->manipFrame.removeFromMouseGrabberPool();
}

void ControlPoint::setState(GrabberState newState)
{
    this->state = newState;
}

void ControlPoint::updateSphere()
{
    if (this->useManipFrame)
        this->position = Vector3(this->manipFrame.position().x, this->manipFrame.position().y, this->manipFrame.position().z);
    this->shape.position = this->position;
    this->shape.buildVerticesFlat();
    this->fromArray(this->shape.mesh.vertexArrayFloat);
    this->update();
}

void ControlPoint::display()
{
    if (this->state == HIDDEN) return;
    else
        this->shader->setVector("color", ControlPoint::GrabberStateColor[this->state]);
    this->updateSphere();
    Mesh::display();
}

void ControlPoint::move(Vector3 newPos)
{
    if (this->useManipFrame) {
        this->manipFrame.setPosition(this->position.x, this->position.y, this->position.z);
    }
    this->position = newPos;
    this->updateSphere();
}
