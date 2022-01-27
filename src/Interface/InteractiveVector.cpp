#include "InteractiveVector.h"

std::shared_ptr<Shader> InteractiveVector::base_shader = nullptr;

InteractiveVector::InteractiveVector()
    : InteractiveVector(Vector3())
{

}
InteractiveVector::InteractiveVector(Vector3 end)
    : InteractiveVector(Vector3(), end)
{

}
InteractiveVector::InteractiveVector(Vector3 start, Vector3 end)
{
    this->startingControlPoint = new ControlPoint(start);
    this->endingControlPoint = new ControlPoint(end);
    this->endingControlPoint->mesh.shader = this->startingControlPoint->mesh.shader;
    this->arrowMesh.shader = this->startingControlPoint->mesh.shader;

    QObject::connect(this->startingControlPoint, &ControlPoint::modified, this, [=](){
        Q_EMIT this->modified(this->getResultingVector());
        Q_EMIT this->startingModified(this->getStartingVector()); });
    QObject::connect(this->endingControlPoint, &ControlPoint::modified, this, [=](){
        Q_EMIT this->modified(this->getResultingVector());
        Q_EMIT this->startingModified(this->getEndingVector()); });
}

void InteractiveVector::display()
{
    this->startingControlPoint->display();
    this->endingControlPoint->display();
    this->arrowMesh.fromArray(this->getArrowPath());
    this->arrowMesh.display(GL_LINES);
}

std::vector<Vector3> InteractiveVector::getArrowPath()
{
    Vector3 start = startingControlPoint->position, end = endingControlPoint->position;
    float length = (end - start).norm();
    float arrowHeadLength = length / 10.f;
    // Todo : make a real arrow with the head
    return {
        start, end,
        end, end + ((start - end).cross(Vector3(0, 0, 1)).normalized() + (start - end).normalized()).normalized() * arrowHeadLength,
        end, end + ((start - end).cross(Vector3(0, 0,-1)).normalized() + (start - end).normalized()).normalized() * arrowHeadLength,
        end, end + ((start - end).cross(Vector3(0, 1, 0)).normalized() + (start - end).normalized()).normalized() * arrowHeadLength,
        end, end + ((start - end).cross(Vector3(0,-1, 0)).normalized() + (start - end).normalized()).normalized() * arrowHeadLength,
    };
}
