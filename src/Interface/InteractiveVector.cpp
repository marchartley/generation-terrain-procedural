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
    this->startingControlPoint = std::make_unique<ControlPoint>(start, 5.f);
    this->endingControlPoint = std::make_unique<ControlPoint>(end, 3.f);
    this->endingControlPoint->mesh.shareShader(this->startingControlPoint->mesh);
    this->startingControlPoint->allowAllAxisRotations(true);
    this->startingControlPoint->allowAllAxisTranslation(true);
    this->endingControlPoint->allowAllAxisTranslation(true);
    this->arrowMesh.shareShader(this->startingControlPoint->mesh);

    QObject::connect(this->startingControlPoint.get(), &ControlPoint::modified, this, [=](){
        Q_EMIT this->modified(this->getResultingVector());
        Q_EMIT this->startingModified(this->getStartingVector());
    });
    QObject::connect(this->endingControlPoint.get(), &ControlPoint::modified, this, [=](){
        Q_EMIT this->modified(this->getResultingVector());
        Q_EMIT this->endingModified(this->getEndingVector());
    });
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
    Vector3 movement = start - this->startingControlPoint->getPosition();
    this->setPositions(start, this->endingControlPoint->getPosition() + movement);
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
    Vector3 start = startingControlPoint->getPosition(), end = endingControlPoint->getPosition();
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
