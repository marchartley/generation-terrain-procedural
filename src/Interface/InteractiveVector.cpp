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
    this->endingControlPoint->mesh.shareShader(this->startingControlPoint->mesh);
    this->arrowMesh.shareShader(this->startingControlPoint->mesh);

    QObject::connect(this->startingControlPoint, &ControlPoint::modified, this, [=](){
        Q_EMIT this->modified(this->getResultingVector());
        Q_EMIT this->startingModified(this->getStartingVector());
    });
    QObject::connect(this->endingControlPoint, &ControlPoint::modified, this, [=](){
        Q_EMIT this->modified(this->getResultingVector());
        Q_EMIT this->startingModified(this->getEndingVector()); });
}

void InteractiveVector::display()
{
    this->startingControlPoint->display();
    this->endingControlPoint->display();
    this->arrowMesh.fromArray(this->getArrowPath());
    this->arrowMesh.display(GL_LINES, 3.f);
}

void InteractiveVector::setPosition(Vector3 start)
{
    Vector3 movement = start - this->startingControlPoint->position;
    this->setPositions(start, this->endingControlPoint->position + movement);
}

void InteractiveVector::setPositions(Vector3 start, Vector3 end)
{
    this->startingControlPoint->move(start);
    this->endingControlPoint->move(end);
}

void InteractiveVector::hide()
{
    this->startingControlPoint->hide();
    this->endingControlPoint->hide();
    this->arrowMesh.hide();
    CustomInteractiveObject::hide();
}

void InteractiveVector::show()
{
    this->startingControlPoint->show();
    this->endingControlPoint->show();
    this->arrowMesh.show();
    CustomInteractiveObject::show();
}

std::vector<Vector3> InteractiveVector::getArrowPath()
{
    Vector3 start = startingControlPoint->position, end = endingControlPoint->position;
    float length = (end - start).norm();
    float arrowHeadLength = length / 10.f;
    return {
        start, end,
        end, end + ((start - end).cross(Vector3(0, 0, 1)).normalized() + (start - end).normalized()).normalized() * arrowHeadLength,
        end, end + ((start - end).cross(Vector3(0, 0,-1)).normalized() + (start - end).normalized()).normalized() * arrowHeadLength,
        end, end + ((start - end).cross(Vector3(0, 1, 0)).normalized() + (start - end).normalized()).normalized() * arrowHeadLength,
        end, end + ((start - end).cross(Vector3(0,-1, 0)).normalized() + (start - end).normalized()).normalized() * arrowHeadLength,
    };
}
