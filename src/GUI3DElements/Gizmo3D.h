#ifndef Gizmo3D_H
#define Gizmo3D_H

#include "DataStructure/Vector3.h"
#include "DataStructure/Quaternion.h"
#include "Graphics/Mesh.h"
#include <QGLViewer/camera.h>
// #include "Interface/CustomInteractiveObject.h"
// #include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/mouseGrabber.h>

enum class GizmoAxis {
    None,
    X,
    Y,
    Z
};

enum class GizmoHandleType {
    None,
    Center,
    TranslateX, TranslateY, TranslateZ,
    RotateX, RotateY, RotateZ
};

enum class DragMode {
    None,
    FreeMove,
    TranslateAxis,
    RotateAxis
};

GizmoAxis axisFromHandle(GizmoHandleType h);
Vector3 axisVector(GizmoAxis axis);

struct GizmoTransform {
    Vector3 position = Vector3(0.f, 0.f, 0.f);
    Quaternion rotation;
};

struct GizmoSettings {
    float radius = 1.f;
    float axisLength = 3.f;
    float ringRadius = 2.f;

    bool visible = true;
    bool displayOnTop = true;
    bool keepConstantScreenSize = true;

    float pickPixelTolerance = 8.f;
    float minRadius = 0.05f;
    float maxRadius = 1e6f;
};

struct GizmoCapabilities {
    bool allowTranslateX = false;
    bool allowTranslateY = false;
    bool allowTranslateZ = false;

    bool allowRotateX = false;
    bool allowRotateY = false;
    bool allowRotateZ = false;

    bool allowFreeMove = true;
};

struct GizmoVisualState {
    GizmoHandleType hovered = GizmoHandleType::None;
    GizmoHandleType active = GizmoHandleType::None;
    Vector3 body_color;

    bool isDragging = false;
};

struct DragState {
    // DragState() = default;

    DragMode mode = DragMode::None;
    GizmoAxis axis = GizmoAxis::None;
    GizmoHandleType handle = GizmoHandleType::None;

    bool active = false;

    Vector3 pressMousePos;

    Vector3 startWorldPosition;
    Quaternion startWorldRotation;

    Vector3 dragPlaneOrigin;
    Vector3 dragPlaneNormal;

    Vector3 startHitPointWorld;
    Vector3 currentHitPointWorld;

    Vector3 startAxisDirectionWorld;
};

struct GizmoState {
    GizmoTransform transform;
    GizmoSettings settings;
    GizmoCapabilities capabilities;
    GizmoVisualState visual;
    DragState drag;
    bool dirty = true;
};

class GizmoInteractor
{
public:
    GizmoInteractor();

    GizmoHandleType pickHandle(const GizmoState& state, const Vector3& mousePos, const qglviewer::Camera* camera) const;

    bool beginDrag(GizmoState& state, const Vector3& mousePos, const qglviewer::Camera* camera);
    bool updateDrag(GizmoState& state, const Vector3& mousePos, const qglviewer::Camera* camera);
    void endDrag(GizmoState& state);

    void updateHover(GizmoState& state, const Vector3& mousePos, const qglviewer::Camera* camera) const;

// private:
    bool computeMouseRay(const Vector3& mousePos, const qglviewer::Camera* camera, Vector3& rayOrigin, Vector3& rayDir) const;

    Vector3 computeTranslationPlaneNormal(GizmoAxis axis, const Vector3& cameraDir) const;

    bool intersectDragPlane(const DragState& drag, const Vector3& rayOrigin, const Vector3& rayDir, Vector3& hitPoint) const;

    float signedAngleAroundAxis(const Vector3& from, const Vector3& to, const Vector3& axis) const;

    bool pickCenter(const GizmoState& state, const Vector3& mousePos, const qglviewer::Camera* camera) const;
    bool pickTranslationAxis(const GizmoState& state, GizmoAxis axis, const Vector3& mousePos, const qglviewer::Camera* camera) const;
    bool pickRotationRing(const GizmoState& state, GizmoAxis axis, const Vector3& mousePos, const qglviewer::Camera* camera) const;
};

class GizmoRenderer
{
public:
    GizmoRenderer();

    void render(GizmoState& state);

private:
    void renderCenter(const GizmoState& state);
    void renderTranslationAxis(const GizmoState& state, GizmoAxis axis);
    void renderRotationRing(const GizmoState& state, GizmoAxis axis);

    std::vector<Vector3> buildCirclePoints(const Vector3& center, GizmoAxis axis, float radius, int segments = 64) const;
    std::vector<Vector3> buildAxisLine(const Vector3& center, GizmoAxis axis, float halfLength) const;
    std::vector<float> colorForHandle(const GizmoState& state, GizmoHandleType handle) const;

    Mesh centerMesh;
    Mesh translationMesh;
    Mesh rotationMesh;
};



class Gizmo3D : /*public QObject,*/ public qglviewer::MouseGrabber {
// Q_OBJECT
public:
    Gizmo3D(bool applyManipulations = true, QObject* parent = nullptr);

    GizmoState& state();
    const GizmoState& state() const;

    void render();

    void mousePressEvent      (QMouseEvent* const event, qglviewer::Camera* const camera);
    void mouseMoveEvent       (QMouseEvent* const event, qglviewer::Camera* const camera);
    void mouseReleaseEvent    (QMouseEvent* const event, qglviewer::Camera* const camera);
    void mouseDoubleClickEvent(QMouseEvent* const event, qglviewer::Camera* const camera);
    void wheelEvent           (QWheelEvent* const event, qglviewer::Camera* const camera);

    bool isManipulated() const;

    void checkIfGrabsMouse(int x, int y, const qglviewer::Camera* const);

    // void setOnPressed(const std::function<void(void)>& onPress) { this->onPointPressedCallbacks.push_back(onPress); }
    // void setOnModified(const std::function<void(void)>& onModified) { this->onPointModifiedCallbacks.push_back(onModified); }
    // void setOnReleased(const std::function<void(void)>& onRelease) { this->onPointReleasedCallbacks.push_back(onRelease); }
    // void setOnTranslation(const std::function<void(void)>& onTranslation) { this->onPointTranslatedCallbacks.push_back(onTranslation); }
    // void setOnRotation(const std::function<void(void)>& onRotation) { this->onPointRotatedCallbacks.push_back(onRotation); }



// Q_SIGNALS:
    // void pointModified();
    // void pointReleased();

    // void translationApplied(const Vector3& translation);
    // void rotationApplied(const Vector3& rotation);

protected:
    GizmoState state_;
    GizmoInteractor interactor_;
    GizmoRenderer renderer_;

    bool applyManipulations;

    DECLARE_EVENT(PointPressed, (), ())
    DECLARE_EVENT(PointModified, (), ())
    DECLARE_EVENT(PointReleased, (), ())
    DECLARE_EVENT(PointTranslated, (), ())
    DECLARE_EVENT(PointRotated, (), ())

    // void emitPressedEvent() { for (auto& fn : onPointPressedCallbacks) { fn(); } }
    // void emitModifiedEvent() { for (auto& fn : onPointModifiedCallbacks) { fn(); } }
    // void emitReleasedEvent() { for (auto& fn : onPointReleasedCallbacks) { fn(); } }
    // void emitTranslatedEvent() { for (auto& fn : onPointTranslatedCallbacks) { fn(); } }
    // void emitRotatedEvent() { for (auto& fn : onPointRotatedCallbacks) { fn(); } }
};

#endif // Gizmo3D_H
