#include "Slider3D.h"

#include "Utils/Utils.h"

Slider3D::Slider3D() : Slider3D(Vector3(), 1.f, 0.f, 0.f, 1.f, Slider3DOrientation::X)
{

}

Slider3D::Slider3D(Vector3 positionMin, float length, float val, float minValue, float maxValue, Slider3DOrientation orientation)
{
    Vector3 positionMax = positionMin + Vector3((orientation & X ? 1 : 0), (orientation & Y ? 1 : 0), (orientation & Z ? 1 : 0)).normalized() * length;
    this->init(positionMin, positionMax, minValue, maxValue, val);
}

Slider3D::Slider3D(Vector3 positionMin, Vector3 positionMax, float val, float minValue, float maxValue)
{
    this->init(positionMin, positionMax, minValue, maxValue, val);
}

void Slider3D::setPosition(Vector3 newPos)
{
    Vector3 movement = newPos - this->minPos;
    this->minPos += movement;
    this->maxPos += movement;
    this->sliderControlPoint->move(this->sliderControlPoint->getPosition() + movement);
    this->sliderControlPoint->custom_constraint = new SliderConstraint(minPos, maxPos);
    this->sliderMesh.fromArray({minPos, maxPos});
}

void Slider3D::setPositions(Vector3 newStart, Vector3 newEnd)
{
    float currentValue = this->getValue();
    this->minPos = newStart;
    this->maxPos = newEnd;
    this->sliderControlPoint->move(Vector3::lerp(currentValue, newStart, newEnd));
    this->constraint = new SliderConstraint(minPos, maxPos);
    this->sliderControlPoint->custom_constraint = constraint;
//    this->sliderControlPoint->setConstraint(new SliderConstraint(newStart, newEnd));
//    std::cout << "New positions : " << newStart << " " << newEnd << std::endl;
    this->sliderMesh.fromArray({minPos, maxPos});
}

void Slider3D::hide()
{
    this->sliderControlPoint->hide();
    this->sliderMesh.hide();
    CustomInteractiveObject::hide();
}

void Slider3D::show()
{
    this->sliderControlPoint->show();
    this->sliderMesh.show();
    CustomInteractiveObject::show();
}

float Slider3D::setValue(float newValue)
{
    // float t = interpolation::linear(newValue, this->minValue, this->maxValue);
    this->sliderControlPoint->move(remap(newValue, minValue, maxValue, minPos, maxPos));
    return newValue;
}

float Slider3D::setValue(Vector3 newPos)
{
    this->sliderControlPoint->move(newPos);
    return this->getValue();
}

void Slider3D::display()
{
    this->sliderControlPoint->display();
    this->sliderMesh.display(GL_LINES, 3.f);
}

float Slider3D::getValue()
{
    return Vector3::remap(this->sliderControlPoint->getPosition(), this->minPos, this->maxPos, this->minValue, this->maxValue);
}

void Slider3D::init(Vector3 positionMin, Vector3 positionMax, float minValue, float maxValue, float val)
{
    this->minPos = positionMin;
    this->maxPos = positionMax;
    this->minValue = minValue;
    this->maxValue = maxValue;
    this->sliderControlPoint = new ControlPoint(remap(val, minValue, maxValue, minPos, maxPos), 5.f);
    this->constraint = new SliderConstraint(positionMin, positionMax);
    this->sliderControlPoint->custom_constraint = constraint;
    this->sliderMesh.fromArray({minPos, maxPos});
    this->sliderMesh.shareShader(this->sliderControlPoint->mesh);

    QObject::connect(this->sliderControlPoint, &ControlPoint::modified, this, [=]() { Q_EMIT this->valueChanged(this->getValue()); });
}


SliderConstraint::SliderConstraint()
{
    this->constraint = new qglviewer::WorldConstraint();
}
SliderConstraint::SliderConstraint(Vector3 minPos, Vector3 maxPos)
    : minPos(minPos), maxPos(maxPos)
{
    this->constraint = new qglviewer::WorldConstraint();
    this->constraint->setTranslationConstraintType(qglviewer::AxisPlaneConstraint::AXIS);
    this->constraint->setRotationConstraintType(qglviewer::AxisPlaneConstraint::FORBIDDEN);
    Vector3 dir = maxPos - minPos;
    this->constraint->setTranslationConstraintDirection(qglviewer::Vec(dir.x, dir.y, dir.z));
}

void SliderConstraint::constrainTranslation(qglviewer::Vec& t, qglviewer::Frame* const fr) {
    this->constraint->constrainTranslation(t, fr);
}

void SliderConstraint::constrainRotation(qglviewer::Quaternion &q, qglviewer::Frame * const fr)
{
    this->constraint->constrainRotation(q, fr);
}
