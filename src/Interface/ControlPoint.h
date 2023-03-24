#ifndef CONTROLPOINT_H
#define CONTROLPOINT_H

#include "Graphics/Sphere.h"
#include "Graphics/Mesh.h"
#include "DataStructure/Vector3.h"

#include "Interface/CustomInteractiveObject.h"
#include <QGLViewer/manipulatedFrame.h>

enum GrabberState {
    HIDDEN         = 0x0,
    INACTIVE       = 0x1,
    ACTIVE         = 0x2,
    POSITIVE       = 0x3,
    NEGATIVE       = 0x4,
    NEUTRAL        = 0x5,
    CUSTOM_STATE_0 = 0x6,
    CUSTOM_STATE_1 = 0x7,
    CUSTOM_STATE_2 = 0x8,
    CUSTOM_STATE_3 = 0x9,
    CUSTOM_STATE_4 = 0xA,
    CUSTOM_STATE_5 = 0xB,
    CUSTOM_STATE_6 = 0xC,
    CUSTOM_STATE_7 = 0xD,
    CUSTOM_STATE_8 = 0xE,
    CUSTOM_STATE_9 = 0xF,
};

/*
class ControlPoint : public qglviewer::ManipulatedFrame
{
    Q_OBJECT
public:
    ControlPoint();
    ControlPoint(Vector3 pos, float radius = 1.f, GrabberState state = INACTIVE, bool useTheManipulatedFrame = true);
    ~ControlPoint();

protected:
    enum ControlPointAction { TRANSLATE, ROTATE, SCALE, NO_ACTION };
    enum Axis {X, Y, Z, NO_AXIS };

public:
    void display();
    void move(Vector3 newPos);

    Vector3 getPosition() const;

    void allowAllAxisTranslation(bool allowed);
    void allowAllAxisRotations(bool allowed);
    void allowAllAxisScaling(bool allowed);

    void mousePressEvent(QMouseEvent* const event  , qglviewer::Camera* const cam );
    void mouseReleaseEvent( QMouseEvent* const event, qglviewer::Camera* const cam);
    void mouseMoveEvent(QMouseEvent* const event, qglviewer::Camera* const cam);
    void wheelEvent(QWheelEvent *const event, qglviewer::Camera *const camera);

    void setGrabberStateColor(std::map<GrabberState, std::vector<float>> stateColorMap);
    void setGrabberStateColor(GrabberState state, std::vector<float> color);

    void setPosition(Vector3 newPosition);
    void setRadius(float newRadius);
    void setState(GrabberState newState);

    Vector3 getFluidTranslation() const;
    Vector3 getFluidRotation() const;
    Vector3 getCurrentTranslation() const;
    Vector3 getCurrentRotation() const;

    qglviewer::Constraint* custom_constraint = nullptr;

    std::map<Axis, bool> allowedTranslations;
    std::map<Axis, bool> allowedRotations;
    std::map<Axis, bool> allowedScaling;

    bool displayOnTop = true;

//    static std::shared_ptr<Shader> base_shader;
    static std::map<GrabberState, std::vector<float>> default_GrabberStateColor;

public Q_SLOTS:
    void hide();
    void show();

Q_SIGNALS:
    void modified();
    void released();

    void translationApplied(Vector3);
    void rotationApplied(Vector3);

protected:
//    void updateStateDependingOnManipFrame();
    void checkIfGrabsMouse(int x, int y,const qglviewer::Camera* const cam);

    Vector3 intersectionWithTranslationWidget(Vector3 rayOrigin, Vector3 rayDir);
    Vector3 intersectionWithRotationWidget(Vector3 rayOrigin, Vector3 rayDir);
    Vector3 intersectionWithScalingWidget(Vector3 rayOrigin, Vector3 rayDir);

    ControlPointAction hoveredAction(Vector3 rayOrigin, Vector3 rayDir);
    Axis hoveredAxis(Vector3 rayOrigin, Vector3 rayDir);

    ControlPointAction currentAction = NO_ACTION;
    Axis currentAxis = NO_AXIS;

    Vector3 position;
    Vector3 mousePositionWhenActionStarted;

    float radius;
    GrabberState currentState = INACTIVE;

//    std::map<GrabberState, std::vector<float>> GrabberStateColor;

//    Mesh translationMesh, rotationMesh, scalingMesh;

    std::map<GrabberState, std::vector<float>> GrabberStateColor;

    std::vector<Vector3> translationHistory;
    std::vector<Vector3> rotationHistory;

    Vector3 getIntersectionWithTranslationAxis(Vector3 rayOrigin, Vector3 rayDir, Axis axis);
    Vector3 getIntersectionWithRotationAxis(Vector3 rayOrigin, Vector3 rayDir, Axis axis);
    Vector3 getIntersectionWithScalingAxis(Vector3 rayOrigin, Vector3 rayDir, Axis axis);

};
*/

class ControlPoint : public qglviewer::ManipulatedFrame
{
    Q_OBJECT
public:
    ControlPoint();
    ControlPoint(Vector3 pos, float radius = 1.f, GrabberState state = INACTIVE, bool useTheManipulatedFrame = true);
    ~ControlPoint();

    void setState(GrabberState newState);
    void updateStateDependingOnManipFrame();

    void checkIfGrabsMouse(int x, int y,const qglviewer::Camera* const cam);

//    void setPosition (const qglviewer::Vec &position);

Q_SIGNALS:
    void modified();
    void released();

    void translationApplied(Vector3);
    void rotationApplied(Vector3);

public Q_SLOTS:
    void hide();
    void show();

private:
    void onUpdate(std::function<void()> func);
    void afterUpdate(std::function<void()> func);
public:
    void updateSphere();
    void move(Vector3 newPos);
    void display();

    void setGrabberStateColor(std::map<GrabberState, std::vector<float>> stateColorMap);
    void setGrabberStateColor(GrabberState state, std::vector<float> color);

    Vector3 getRotation() const { return Vector3::quaternionToEuler(this->rotation()); }
    Vector3 getPosition() const { return this->position(); }
    Vector3 getFluidTranslation() const {
        if (positionsHistory.empty()) return Vector3();
        return (this->getPosition() - this->positionsHistory.front()).normalize(); };
    Vector3 getLastMovement() const { return (this->prevPosition - this->getPosition()).normalize(); };

    void mousePressEvent(QMouseEvent* const event  , qglviewer::Camera* const cam );
    void mouseReleaseEvent( QMouseEvent* const event, qglviewer::Camera* const cam);
    void mouseMoveEvent(QMouseEvent* const event, qglviewer::Camera* const cam);
    void wheelEvent(QWheelEvent *const event, qglviewer::Camera *const camera);

    void setSphereRadius(float newRadius);

//    Vector3 position;
//    Vector3 pos;
    std::vector<Vector3> positionsHistory;
    Vector3 prevPosition;
    GrabberState state;
    bool useManipFrame;
    bool currentlyManipulated;

    Vector3 initialPosition;
    Vector3 initialRotation;

    Vector3 pressedPosBeforeAction;

    Mesh mesh;
    Sphere shape;

    Mesh translationMeshes;
    Mesh rotationMeshes;
    Mesh rotationHelperSphere;

    float arrowSize;
    float circleRadius;

    void setRadius(float newRadius) { this->radius = newRadius; }
    float radius;
    float minSphereRadius = -1;
    float maxSphereRadius = -1;

    std::function<void()> onUpdateCallback;
    std::function<void()> afterUpdateCallback;

    qglviewer::ManipulatedFrame manipFrame;
    std::map<GrabberState, std::vector<float>> GrabberStateColor;

    static std::shared_ptr<Shader> base_shader;
    static std::map<GrabberState, std::vector<float>> default_GrabberStateColor;

    void allowAllAxisTranslation(bool allow);
    void allowAllAxisRotations(bool allow);

    enum Axis {X, Y, Z, NONE};
    std::map<Axis, bool> allowedTranslations;
    std::map<Axis, bool> allowedRotations;

    bool isApplyingFreeMove = false;
    bool isApplyingTranslation = false;
    bool isApplyingRotation = false;
    Axis currentAxis = NONE;

    bool displayOnTop = true;

    qglviewer::WorldConstraint* default_constraint;
    qglviewer::Constraint* custom_constraint = nullptr;
protected:
    std::vector<Vector3> computeCircle(Axis axis);
    Vector3 getIntersectionWithPlane(Vector3 rayOrigin, Vector3 rayDir, Axis axis);

    bool mouseOnCentralSphere(Vector3 rayOrigin, Vector3 rayDir);
    bool mouseOnTranslationArrow(Vector3 rayOrigin, Vector3 rayDir);
    bool mouseOnRotationCircle(Vector3 rayOrigin, Vector3 rayDir);

    bool stillOnInitialState = true;

};


class CustomConstraint : public qglviewer::Constraint
{
public:
    CustomConstraint();

    virtual void constrainTranslation(qglviewer::Vec& t, qglviewer::Frame* const fr);
    virtual void constrainRotation(qglviewer::Quaternion& q, qglviewer::Frame* const fr);

private:
    qglviewer::WorldConstraint* constraint;

    bool useTranslation;
};


#endif // CONTROLPOINT_H
