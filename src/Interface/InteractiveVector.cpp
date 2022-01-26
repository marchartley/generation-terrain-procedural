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
    this->startingControlPoint = ControlPoint(start);
    this->endingControlPoint = ControlPoint(end);
}

void InteractiveVector::display()
{
    this->arrowMesh.fromArray(this->getArrowPath());
    this->arrowMesh.display(GL_LINES);
    this->startingControlPoint.display();
    this->endingControlPoint.display();
}

void InteractiveVector::onUpdate(std::function<void ()> callback)
{
    this->onUpdateCallback = callback;

    this->startingControlPoint.afterUpdate(callback);
    this->endingControlPoint.afterUpdate(callback);
}

std::vector<Vector3> InteractiveVector::getArrowPath()
{
    Vector3 start = startingControlPoint.position, end = endingControlPoint.position;
//    float length = (end - start).norm();
//    float arrowHeadLength = length / 10.f;
    // Todo : make a real arrow with the head
    return {
        start, end
    };
}
