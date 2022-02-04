#ifndef SLIDER3D_H
#define SLIDER3D_H

#include "Interface/ControlPoint.h"
#include "Interface/CustomInteractiveObject.h"
#include "DataStructure/Vector3.h"

enum Slider3DOrientation {
    X = 0b001,
    Y = 0b010,
    Z = 0b100
};

class SliderConstraint;
class Slider3D : public CustomInteractiveObject
{
    Q_OBJECT
public:
    Slider3D();
    Slider3D(Vector3 positionMin, float length, float val = 0.f, float minValue = 0.f, float maxValue = 1.f, Slider3DOrientation orientation = X);
    Slider3D(Vector3 positionMin, Vector3 positionMax, float val = 0.f, float minValue = 0.f, float maxValue = 1.f);

public Q_SLOTS:
    void hide();
    void show();

Q_SIGNALS:
    void valueChanged(float newValue);

public:
    void display();

    float getValue();

    float minValue;
    float maxValue;
    Vector3 minPos;
    Vector3 maxPos;

    ControlPoint *sliderControlPoint;
    Mesh sliderMesh;

    void init(Vector3 positionMin, Vector3 positionMax, float minValue = 0.f, float maxValue = 1.f, float val = 0.f);
};


class SliderConstraint : public qglviewer::Constraint
{
public:
    SliderConstraint();
    SliderConstraint(Vector3 minPos, Vector3 maxPos);

    virtual void constrainTranslation(qglviewer::Vec& t, qglviewer::Frame* const fr);
    virtual void constrainRotation(qglviewer::Quaternion& q, qglviewer::Frame* const fr);

private:
    qglviewer::AxisPlaneConstraint* constraint;

    Vector3 minPos;
    Vector3 maxPos;
};

#endif // SLIDER3D_H
