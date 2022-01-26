#include "Slider3D.h"

#include "Utils/Utils.h"

Slider3D::Slider3D()
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

void Slider3D::display()
{
    this->sliderMesh.display(GL_LINES);
    this->sliderControlPoint.display();
}

float Slider3D::getValue()
{
    return Vector3::remap(this->sliderControlPoint.position, this->minPos, this->maxPos, this->minValue, this->maxValue);
}

void Slider3D::init(Vector3 positionMin, Vector3 positionMax, float minValue, float maxValue, float val)
{
    this->minPos = positionMin;
    this->maxPos = positionMax;
    this->minValue = minValue;
    this->maxValue = maxValue;
    this->sliderControlPoint = ControlPoint(remap(val, minValue, maxValue, minPos, maxPos), 5.f);
    this->sliderControlPoint.manipFrame.setConstraint(new SliderConstraint(positionMin, positionMax));
    this->sliderMesh.fromArray({minPos, maxPos});
}


SliderConstraint::SliderConstraint()
{
    this->constraint = new qglviewer::LocalConstraint();
}
SliderConstraint::SliderConstraint(Vector3 minPos, Vector3 maxPos)
    : minPos(minPos), maxPos(maxPos)
{
    this->constraint = new qglviewer::LocalConstraint();
    this->constraint->setTranslationConstraintType(qglviewer::AxisPlaneConstraint::AXIS);
    Vector3 dir = maxPos - minPos;
    this->constraint->setTranslationConstraintDirection(qglviewer::Vec(dir.x, dir.y, dir.z));
}

void SliderConstraint::constrainTranslation(qglviewer::Vec& t, qglviewer::Frame* const fr) {
    this->constraint->constrainTranslation(t, fr);/*
    if (std::abs(t.x) > std::abs(t.y) && std::abs(t.x) > std::abs(t.z)) {
        this->constraintx->constrainTranslation(t, fr);
    } else if (std::abs(t.y) > std::abs(t.x) && std::abs(t.y) > std::abs(t.z)) {
        this->constrainty->constrainTranslation(t, fr);
    } else if (std::abs(t.z) > std::abs(t.x) && std::abs(t.z) > std::abs(t.y)) {
        this->constraintz->constrainTranslation(t, fr);
    }*/
}
