#ifndef CONTROLPOINT_H
#define CONTROLPOINT_H

// #define USE_OLD_CONTROL_POINTS 0

// #if !USE_OLD_CONTROL_POINTS


#include "Graphics/Sphere.h"
#include "Graphics/Mesh.h"
#include "DataStructure/Vector3.h"

#include "Interface/CustomInteractiveObject.h"
#include <QGLViewer/manipulatedFrame.h>
#include "GUI3DElements/Gizmo3D.h"

class ControlPoint : public Gizmo3D {
// Q_OBJECT
public:
    enum GrabberState {
        HIDDEN, DEFAULT, POSITIVE, NEGATIVE, NEUTRAL, CUSTOM_STATE_0, CUSTOM_STATE_1, CUSTOM_STATE_2, CUSTOM_STATE_3, CUSTOM_STATE_4, CUSTOM_STATE_5, CUSTOM_STATE_6, CUSTOM_STATE_7, CUSTOM_STATE_8, CUSTOM_STATE_9
    };

    ControlPoint();
    ControlPoint(const Vector3& pos, float radius = 1.f, GrabberState state = DEFAULT, bool applyManipulations = true);

    Vector3 getPosition() const;
    Vector3 getRotation() const;

    void allowAllAxisRotations(bool allow = true);
    void allowAllAxisTranslation(bool allow = true);

    void display();

    void hide();
    void show();
    void setVisible(bool visibility);

    void move(const Vector3& newPos);
    inline void setPosition(const Vector3& newPos) { move(newPos); }

    void setGrabberStateColor(GrabberState state, Vector3 color);
    void setGrabberStateColor(std::map<GrabberState, Vector3> stateColorMap);

    void setState(GrabberState newState);

    void setRadius(float newRadius);

    void setConstraint(qglviewer::Constraint* constraint);

    void setDisplayOnTop(bool enable);

    Vector3 getFluidTranslation();

protected:
    std::map<GrabberState, Vector3> stateColorMap;


protected:
    static std::map<GrabberState, Vector3> default_GrabberStateColor;
};




/*
struct InteractionState;
struct RenderState;
struct ConstraintState;


enum class ControlPointAxis {X, Y, Z, NONE};

struct InteractionState {
    bool isApplyingFreeMove = false;
    bool isApplyingTranslation = false;
    bool isApplyingRotation = false;
    ControlPointAxis currentAxis = ControlPointAxis::NONE;
    Vector3 pressedPosBeforeAction;
    Vector3 currentMousePosOnAction;
};

class ControlPoint__ : public qglviewer::ManipulatedFrame
{
    Q_OBJECT
public:

    enum GrabberState {
        HIDDEN, INACTIVE, ACTIVE, POSITIVE, NEGATIVE, NEUTRAL, CUSTOM_STATE_0, CUSTOM_STATE_1, CUSTOM_STATE_2, CUSTOM_STATE_3, CUSTOM_STATE_4, CUSTOM_STATE_5, CUSTOM_STATE_6, CUSTOM_STATE_7, CUSTOM_STATE_8, CUSTOM_STATE_9
    };

    ControlPoint__();
    ControlPoint__(const Vector3& pos, float radius = 1.f, GrabberState state = INACTIVE, bool useTheManipulatedFrame = true);
    ~ControlPoint__();

    void setState(GrabberState newState);
    void setVisible(bool visibility);
    void updateStateDependingOnManipFrame();
    void setDisplayOnTop(bool enabled) { this->displayOnTop = enabled; }

    void checkIfGrabsMouse(int x, int y,const qglviewer::Camera* const cam);

    void updateSphere();
    void move(const Vector3& newPos);
    void display();

    void setGrabberStateColor(std::map<GrabberState, std::vector<float>> stateColorMap);
    void setGrabberStateColor(GrabberState state, std::vector<float> color);

    void setPosition(const Vector3& newPos);
    void setPosition(float x, float y, float z);

    Vector3 getRotation() const; // { return Vector3::quaternionToEuler(this->rotation()); }
    Vector3 getPosition() const;
    Vector3 getFluidTranslation() const;
    Vector3 getLastMovement() const;

    void mousePressEvent(QMouseEvent* const event  , qglviewer::Camera* const cam );
    void mouseReleaseEvent( QMouseEvent* const event, qglviewer::Camera* const cam);
    void mouseMoveEvent(QMouseEvent* const event, qglviewer::Camera* const cam);
    void wheelEvent(QWheelEvent *const event, qglviewer::Camera *const camera);

    void setSphereRadius(float newRadius);

    void setRadius(float newRadius);

    void allowAllAxisTranslation(bool allow);
    void allowAllAxisRotations(bool allow);

    void setConstraint(qglviewer::Constraint * const constraint);




    bool displayOnTop = true;

Q_SIGNALS:
    void pointModified();
    void pointReleased();

    void translationApplied(const Vector3&);
    void rotationApplied(const Vector3&);

public Q_SLOTS:
    void hide();
    void show();

protected:
    void onUpdate(std::function<void()> func);
    void afterUpdate(std::function<void()> func);

    void updateInteractionState();
    void updateGeometryIfNeeded();

//    Vector3 position;
//    Vector3 pos;
    std::vector<Vector3> positionsHistory;
    Vector3 prevPosition;
    Vector3 currentPosition;
    GrabberState state;
    bool useManipFrame = false;
    bool currentlyManipulated = false;

    Vector3 initialPosition;
    Vector3 initialRotation;
    Mesh mesh;

    Mesh translationMeshes;
    Mesh rotationMeshes;
    Mesh rotationHelperSphere;

    float arrowSize = 1.f;
    float circleRadius = 1.f;

    float radius = 1.f;
    float minSphereRadius = -1;
    float maxSphereRadius = -1;

    std::function<void()> onUpdateCallback;
    std::function<void()> afterUpdateCallback;

    std::map<GrabberState, std::vector<float>> GrabberStateColor;

    static std::shared_ptr<Shader> base_shader;
    static std::map<GrabberState, std::vector<float>> default_GrabberStateColor;

    std::map<ControlPointAxis, bool> allowedTranslations;
    std::map<ControlPointAxis, bool> allowedRotations;

    InteractionState interactionState;

    qglviewer::WorldConstraint* default_constraint = nullptr;
    qglviewer::Constraint* custom_constraint = nullptr;

    std::vector<Vector3> computeCircle(ControlPointAxis axis);
    Vector3 getIntersectionWithPlane(const Vector3& rayOrigin, const Vector3& rayDir, ControlPointAxis axis);

    bool mouseOnCentralSphere(const Vector3& rayOrigin, const Vector3& rayDir);
    bool mouseOnTranslationArrow(const Vector3& rayOrigin, const Vector3& rayDir);
    std::pair<ControlPointAxis, Vector3> mouseOnRotationCircle(const Vector3& rayOrigin, const Vector3& rayDir);

    bool geometryDirty = true;

    bool visible = true;
};



struct RenderState {
    float radius = 1.f;
    float arrowSize = 1.f;
    float circleRadius = 1.f;
    bool geometryDirty = true;
};

struct ConstraintState {
    std::unique_ptr<qglviewer::WorldConstraint> hoverConstraint;
    qglviewer::Constraint* customConstraint = nullptr;
};
*/


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
