#include "Gizmo3D.h"

#include "Utils/Collisions.h"
#include "Graphics/Sphere.h"


std::map<Gizmo3D::GrabberState, std::vector<float>> Gizmo3D::default_GrabberStateColor = {
    {Gizmo3D::GrabberState::HIDDEN, {.0f, .0f, .0f, 0.f}},
    {Gizmo3D::GrabberState::INACTIVE, {.3f, .0f, .0f, .5f}},
    {Gizmo3D::GrabberState::ACTIVE, {.8f, .0f, .0f, .8f}},
    {Gizmo3D::GrabberState::POSITIVE, {.2f, 1.f, .1f, .8f}},
    {Gizmo3D::GrabberState::NEGATIVE, {1.f, .2f, .1f, .8f}},
    {Gizmo3D::GrabberState::NEUTRAL, {.8f, .8f, .8f, .8f}},
    };

GizmoAxis axisFromHandle(GizmoHandleType h) {
    switch (h) {
    case GizmoHandleType::TranslateX:
    case GizmoHandleType::RotateX: return GizmoAxis::X;
    case GizmoHandleType::TranslateY:
    case GizmoHandleType::RotateY: return GizmoAxis::Y;
    case GizmoHandleType::TranslateZ:
    case GizmoHandleType::RotateZ: return GizmoAxis::Z;
    default: return GizmoAxis::None;
    }
}

Vector3 axisVector(GizmoAxis axis) {
    switch (axis) {
    case GizmoAxis::X: return Vector3(1.f, 0.f, 0.f);
    case GizmoAxis::Y: return Vector3(0.f, 1.f, 0.f);
    case GizmoAxis::Z: return Vector3(0.f, 0.f, 1.f);
    default: return Vector3(0.f, 0.f, 0.f);
    }
}

GizmoInteractor::GizmoInteractor()
{

}

GizmoHandleType GizmoInteractor::pickHandle(const GizmoState &state, const Vector3 &mousePos, const qglviewer::Camera *camera) const
{
    if (!state.settings.visible || !camera)
        return GizmoHandleType::None;

    if (state.capabilities.allowFreeMove && pickCenter(state, mousePos, camera))
        return GizmoHandleType::Center;

    if (state.capabilities.allowTranslateX && pickTranslationAxis(state, GizmoAxis::X, mousePos, camera))
        return GizmoHandleType::TranslateX;
    if (state.capabilities.allowTranslateY && pickTranslationAxis(state, GizmoAxis::Y, mousePos, camera))
        return GizmoHandleType::TranslateY;
    if (state.capabilities.allowTranslateZ && pickTranslationAxis(state, GizmoAxis::Z, mousePos, camera))
        return GizmoHandleType::TranslateZ;

    if (state.capabilities.allowRotateX && pickRotationRing(state, GizmoAxis::X, mousePos, camera))
        return GizmoHandleType::RotateX;
    if (state.capabilities.allowRotateY && pickRotationRing(state, GizmoAxis::Y, mousePos, camera))
        return GizmoHandleType::RotateY;
    if (state.capabilities.allowRotateZ && pickRotationRing(state, GizmoAxis::Z, mousePos, camera))
        return GizmoHandleType::RotateZ;

    return GizmoHandleType::None;
}

bool GizmoInteractor::beginDrag(GizmoState &state, const Vector3 &mousePos, const qglviewer::Camera *camera)
{
    if (!camera)
        return false;

    Vector3 rayOrigin, rayDir;
    if (!computeMouseRay(mousePos, camera, rayOrigin, rayDir))
        return false;

    GizmoHandleType handle = this->pickHandle(state, mousePos, camera);

    if (handle == GizmoHandleType::None)
        return false;

    DragState drag;
    drag.active = true;
    drag.handle = handle;
    drag.axis = axisFromHandle(handle);
    drag.pressMousePos = mousePos;
    drag.startWorldPosition = state.transform.position;
    drag.startWorldRotation = state.transform.rotation;

    Vector3 hit;

    if (handle == GizmoHandleType::Center) {
        drag.mode = DragMode::FreeMove;

        // Plane perpendicular to camera, through current position
        drag.dragPlaneOrigin = state.transform.position;
        drag.dragPlaneNormal = rayDir.normalized();

        if (!intersectDragPlane(drag, rayOrigin, rayDir, hit))
            hit = state.transform.position;

        drag.startHitPointWorld = hit;
        drag.currentHitPointWorld = hit;
    }
    else if (handle == GizmoHandleType::TranslateX ||
             handle == GizmoHandleType::TranslateY ||
             handle == GizmoHandleType::TranslateZ) {
        drag.mode = DragMode::TranslateAxis;

        Vector3 camDir = rayDir.normalized();
        drag.dragPlaneOrigin = state.transform.position;
        drag.dragPlaneNormal = computeTranslationPlaneNormal(drag.axis, camDir);

        if (!intersectDragPlane(drag, rayOrigin, rayDir, hit))
            return false;

        drag.startHitPointWorld = hit;
        drag.currentHitPointWorld = hit;
    }
    else if (handle == GizmoHandleType::RotateX ||
             handle == GizmoHandleType::RotateY ||
             handle == GizmoHandleType::RotateZ) {
        drag.mode = DragMode::RotateAxis;

        Vector3 axis = axisVector(drag.axis);
        drag.dragPlaneOrigin = state.transform.position;
        drag.dragPlaneNormal = axis;

        if (!intersectDragPlane(drag, rayOrigin, rayDir, hit))
            return false;

        drag.startHitPointWorld = hit;
        drag.currentHitPointWorld = hit;
    }
    else {
        return false;
    }

    state.drag = drag;
    state.visual.active = handle;
    state.visual.isDragging = true;
    return true;
}

bool GizmoInteractor::updateDrag(GizmoState &state, const Vector3 &mousePos, const qglviewer::Camera *camera)
{
    if (!camera || !state.drag.active)
        return false;

    Vector3 rayOrigin, rayDir;
    if (!computeMouseRay(mousePos, camera, rayOrigin, rayDir))
        return false;

    Vector3 hit;
    if (!intersectDragPlane(state.drag, rayOrigin, rayDir, hit))
        return false;

    state.drag.currentHitPointWorld = hit;

    switch (state.drag.mode) {
    case DragMode::FreeMove: {
        Vector3 delta = state.drag.currentHitPointWorld - state.drag.startHitPointWorld;
        state.transform.position = state.drag.startWorldPosition + delta;
        return true;
    }

    case DragMode::TranslateAxis: {
        Vector3 axis = axisVector(state.drag.axis).normalized();
        Vector3 delta = state.drag.currentHitPointWorld - state.drag.startHitPointWorld;
        float projectedAmount = delta.dot(axis);
        state.transform.position = state.drag.startWorldPosition + axis * projectedAmount;
        return true;
    }

    case DragMode::RotateAxis: {
        Vector3 center = state.drag.dragPlaneOrigin;
        Vector3 axis = axisVector(state.drag.axis).normalized();

        Vector3 v0 = (state.drag.startHitPointWorld - center);
        Vector3 v1 = (state.drag.currentHitPointWorld - center);

        if (v0.norm2() < 1e-8f || v1.norm2() < 1e-8f)
            return false;

        float angle = signedAngleAroundAxis(v0, v1, axis);
        Quaternion dq = Quaternion::AxisAngle(axis, angle);
        state.transform.rotation = dq * state.drag.startWorldRotation;
        return true;
    }

    default:
        return false;
    }
}

void GizmoInteractor::endDrag(GizmoState &state)
{
    state.drag = DragState{};
    state.visual.active = GizmoHandleType::None;
    state.visual.isDragging = false;
}

void GizmoInteractor::updateHover(GizmoState &state, const Vector3 &mousePos, const qglviewer::Camera *camera) const
{
    if (state.drag.active)
        return;

    state.visual.hovered = pickHandle(state, mousePos, camera);
}

bool GizmoInteractor::computeMouseRay(const Vector3 &mousePos, const qglviewer::Camera *camera, Vector3 &rayOrigin, Vector3 &rayDir) const
{
    if (!camera) return false;

    qglviewer::Vec orig, dir;
    camera->convertClickToLine(QPoint(mousePos.x(), mousePos.y()), orig, dir);

    rayOrigin = Vector3(orig);
    rayDir = Vector3(dir).normalized();
    return rayDir.norm2() > 0.f;
}

Vector3 GizmoInteractor::computeTranslationPlaneNormal(GizmoAxis axis, const Vector3 &cameraDir) const
{
    Vector3 a = axisVector(axis).normalized();
    Vector3 c = cameraDir.normalized();

    Vector3 n = a.cross(c.cross(a));
    if (n.norm2() < 1e-6f) {
        // camera almost parallel to axis, fallback
        Vector3 fallback = (std::abs(a.x()) < 0.9f) ? Vector3(1.f, 0.f, 0.f)
                                                    : Vector3(0.f, 1.f, 0.f);
        n = a.cross(fallback.cross(a));
    }

    return n.normalized();
}

bool GizmoInteractor::intersectDragPlane(const DragState &drag, const Vector3 &rayOrigin, const Vector3 &rayDir, Vector3 &hitPoint) const
{
    const Vector3& planeOrigin = drag.dragPlaneOrigin;
    const Vector3& planeNormal = drag.dragPlaneNormal;

    hitPoint = Collision::intersectionRayPlane(rayOrigin, rayDir, planeOrigin, planeNormal);
    return hitPoint.isValid();
}

float GizmoInteractor::signedAngleAroundAxis(const Vector3 &from, const Vector3 &to, const Vector3 &axis) const
{
    return from.getSignedAngleAroundAxisWith(to, axis);
}

bool GizmoInteractor::pickCenter(const GizmoState &state, const Vector3 &mousePos, const qglviewer::Camera *camera) const
{
    Vector3 rayOrigin, rayDir;
    if (!computeMouseRay(mousePos, camera, rayOrigin, rayDir))
        return false;

    return Collision::intersectionRaySphere(rayOrigin, rayDir, state.transform.position, state.settings.radius).isValid();
}

bool GizmoInteractor::pickTranslationAxis(const GizmoState &state, GizmoAxis axis, const Vector3 &mousePos, const qglviewer::Camera *camera) const
{
    if (!camera || axis == GizmoAxis::None)
        return false;

    Vector3 dir = axisVector(axis);
    Vector3 center = state.transform.position;
    float halfLength = state.settings.axisLength;

    Vector3 p0 = center - dir * halfLength;
    Vector3 p1 = center + dir * halfLength;

    qglviewer::Vec s0 = camera->projectedCoordinatesOf(qglviewer::Vec(p0.x(), p0.y(), p0.z()));
    qglviewer::Vec s1 = camera->projectedCoordinatesOf(qglviewer::Vec(p1.x(), p1.y(), p1.z()));

    float d = Collision::shortestDistanceSqrToSegment(mousePos, s0, s1);
    return d <= state.settings.pickPixelTolerance * state.settings.pickPixelTolerance;
}

bool GizmoInteractor::pickRotationRing(const GizmoState &state, GizmoAxis axis, const Vector3 &mousePos, const qglviewer::Camera *camera) const
{
    if (!camera || axis == GizmoAxis::None)
        return false;

    const int segments = 64;
    const Vector3 center = state.transform.position;
    const float radius = state.settings.ringRadius;

    std::vector<Vector3> pts;
    pts.reserve(segments + 1);

    for (int i = 0; i <= segments; ++i) {
        float t = (2.f * float(M_PI) * i) / segments;
        Vector3 p;

        switch (axis) {
        case GizmoAxis::X:
            p = center + Vector3(0.f, std::cos(t), std::sin(t)) * radius;
            break;
        case GizmoAxis::Y:
            p = center + Vector3(std::cos(t), 0.f, std::sin(t)) * radius;
            break;
        case GizmoAxis::Z:
            p = center + Vector3(std::cos(t), std::sin(t), 0.f) * radius;
            break;
        default:
            return false;
        }

        pts.push_back(p);
    }

    float best = std::numeric_limits<float>::max();

    for (int i = 0; i < segments; ++i) {
        qglviewer::Vec s0 = camera->projectedCoordinatesOf(qglviewer::Vec(pts[i].x(), pts[i].y(), pts[i].z()));
        qglviewer::Vec s1 = camera->projectedCoordinatesOf(qglviewer::Vec(pts[i + 1].x(), pts[i + 1].y(), pts[i + 1].z()));

        float d = Collision::shortestDistanceSqrToSegment(mousePos, s0, s1);
        best = std::min(best, d);
    }

    return best <= state.settings.pickPixelTolerance * state.settings.pickPixelTolerance;
}






GizmoRenderer::GizmoRenderer()
{
    this->centerMesh.update();
    this->rotationMesh.update();
    this->translationMesh.update();
}

void GizmoRenderer::render(GizmoState &state)
{
    if (!state.settings.visible)
        return;


    GLboolean m_origin_blend, m_origin_depth, m_origin_cull;
    glGetBooleanv(GL_BLEND, &m_origin_blend);
    glGetBooleanv(GL_DEPTH_TEST,&m_origin_depth);
    glGetBooleanv(GL_CULL_FACE, &m_origin_cull);
    if (state.settings.displayOnTop) {
        glEnable(GL_BLEND);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);
    }
    renderCenter(state);

    if (state.capabilities.allowTranslateX) renderTranslationAxis(state, GizmoAxis::X);
    if (state.capabilities.allowTranslateY) renderTranslationAxis(state, GizmoAxis::Y);
    if (state.capabilities.allowTranslateZ) renderTranslationAxis(state, GizmoAxis::Z);

    if (state.capabilities.allowRotateX) renderRotationRing(state, GizmoAxis::X);
    if (state.capabilities.allowRotateY) renderRotationRing(state, GizmoAxis::Y);
    if (state.capabilities.allowRotateZ) renderRotationRing(state, GizmoAxis::Z);

    if (m_origin_blend == GL_TRUE) glEnable(GL_BLEND);
    else glDisable(GL_BLEND);
    if (m_origin_depth == GL_TRUE) glEnable(GL_DEPTH_TEST);
    else glDisable(GL_DEPTH_TEST);
    if (m_origin_cull == GL_TRUE) glEnable(GL_CULL_FACE);
    else glDisable(GL_CULL_FACE);

    state.dirty = false;
}

void GizmoRenderer::renderCenter(const GizmoState &state)
{
    if (state.dirty) {
        Sphere s(state.settings.radius, state.transform.position, 10, 10);
        centerMesh.fromArray(s.mesh.vertexArrayFloat);
    }
    std::vector<float> color = (state.visual.hovered == GizmoHandleType::Center ? std::vector<float>{1.f, 0.f, 0.f, 1.f} : std::vector<float>{.8f, 0.f, 0.f, 1.f});
    centerMesh.shader->setVector("color", color);
    centerMesh.display();
}

void GizmoRenderer::renderTranslationAxis(const GizmoState &state, GizmoAxis axis)
{
    if (state.dirty) {
        translationMesh.fromArray(buildAxisLine(state.transform.position, axis, state.settings.axisLength));
    }
    std::vector<float> color = (axis == GizmoAxis::X ? std::vector<float>{1, 0, 0, 1} : axis == GizmoAxis::Y ? std::vector<float>{0, 1, 0, 1} : std::vector<float>{0, 0, 1, 1});
    // std::vector<float> color = colorForHandle(state, (axis == GizmoAxis::X ? GizmoHandleType::TranslateX : axis == GizmoAxis::Y ? GizmoHandleType::TranslateY : GizmoHandleType::TranslateZ));
    bool hovered = (axis == GizmoAxis::X && state.visual.hovered == GizmoHandleType::TranslateX ? true : axis == GizmoAxis::Y && state.visual.hovered == GizmoHandleType::TranslateY ? true : axis == GizmoAxis::Z && state.visual.hovered == GizmoHandleType::TranslateZ ? true : false);
    bool active = (axis == GizmoAxis::X && state.visual.active == GizmoHandleType::TranslateX ? true : axis == GizmoAxis::Y && state.visual.active == GizmoHandleType::TranslateY ? true : axis == GizmoAxis::Z && state.visual.active == GizmoHandleType::TranslateZ ? true : false);
    if (!active)
        color[3] = .8f;

    translationMesh.shader->setVector("color", color);
    translationMesh.display(GL_LINES, (hovered ? 5 : 1));
}

void GizmoRenderer::renderRotationRing(const GizmoState &state, GizmoAxis axis)
{
    if (state.dirty) {
        rotationMesh.fromArray(buildCirclePoints(state.transform.position, axis, state.settings.ringRadius));
    }
    std::vector<float> color = (axis == GizmoAxis::X ? std::vector<float>{1, 0, 0, 1} : axis == GizmoAxis::Y ? std::vector<float>{0, 1, 0, 1} : std::vector<float>{0, 0, 1, 1});
    bool hovered = (axis == GizmoAxis::X && state.visual.hovered == GizmoHandleType::RotateX ? true : axis == GizmoAxis::Y && state.visual.hovered == GizmoHandleType::RotateY ? true : axis == GizmoAxis::Z && state.visual.hovered == GizmoHandleType::RotateZ ? true : false);
    bool active = (axis == GizmoAxis::X && state.visual.active == GizmoHandleType::RotateX ? true : axis == GizmoAxis::Y && state.visual.active == GizmoHandleType::RotateY ? true : axis == GizmoAxis::Z && state.visual.active == GizmoHandleType::RotateZ ? true : false);
    if (!active)
        color[3] = .8f;

    rotationMesh.shader->setVector("color", color);
    rotationMesh.display(GL_LINES, (hovered ? 5 : 1));
}

std::vector<Vector3> GizmoRenderer::buildCirclePoints(const Vector3 &center, GizmoAxis axis, float radius, int segments) const
{
    std::vector<Vector3> pts;
    pts.reserve(segments * 2);

    for (int i = 0; i < segments; ++i) {
        float t0 = (2.f * float(M_PI) * i) / segments;
        float t1 = (2.f * float(M_PI) * (i + 1)) / segments;

        Vector3 p0, p1;

        switch (axis) {
        case GizmoAxis::X:
            p0 = center + Vector3(0.f, std::cos(t0), std::sin(t0)) * radius;
            p1 = center + Vector3(0.f, std::cos(t1), std::sin(t1)) * radius;
            break;
        case GizmoAxis::Y:
            p0 = center + Vector3(std::cos(t0), 0.f, std::sin(t0)) * radius;
            p1 = center + Vector3(std::cos(t1), 0.f, std::sin(t1)) * radius;
            break;
        case GizmoAxis::Z:
            p0 = center + Vector3(std::cos(t0), std::sin(t0), 0.f) * radius;
            p1 = center + Vector3(std::cos(t1), std::sin(t1), 0.f) * radius;
            break;
        default:
            continue;
        }

        pts.push_back(p0);
        pts.push_back(p1);
    }

    return pts;
}

std::vector<Vector3> GizmoRenderer::buildAxisLine(const Vector3 &center, GizmoAxis axis, float halfLength) const
{
    Vector3 dir = axisVector(axis);
    return {
        center - dir * halfLength,
        center + dir * halfLength
    };
}

std::vector<float> GizmoRenderer::colorForHandle(const GizmoState &state, GizmoHandleType handle) const
{
    bool active = (state.visual.active == handle);
    bool hovered = (state.visual.hovered == handle);

    if (active) {
        return {1.f, 1.f, 0.f, 1.f}; // yellow
    }
    if (hovered) {
        return {1.f, 1.f, 1.f, 1.f}; // white
    }

    switch (handle) {
    case GizmoHandleType::TranslateX:
    case GizmoHandleType::RotateX:
        return {1.f, 0.f, 0.f, 1.f};

    case GizmoHandleType::TranslateY:
    case GizmoHandleType::RotateY:
        return {0.f, 1.f, 0.f, 1.f};

    case GizmoHandleType::TranslateZ:
    case GizmoHandleType::RotateZ:
        return {0.f, 0.f, 1.f, 1.f};

    case GizmoHandleType::Center:
        return {0.8f, 0.8f, 0.8f, 1.f};

    default:
        return {0.5f, 0.5f, 0.5f, 1.f};
    }
}






Gizmo3D::Gizmo3D(bool applyManipulations, QObject* parent)
    : QObject(parent), applyManipulations(applyManipulations)
{
    if (applyManipulations) {
        this->addInMouseGrabberPool();
    }
}

GizmoState &Gizmo3D::state()
{
    this->state_.dirty = true;
    return this->state_;
}

const GizmoState &Gizmo3D::state() const
{
    return this->state_;
}

void Gizmo3D::render()
{
    this->renderer_.render(state());
}


void Gizmo3D::checkIfGrabsMouse(int x, int y, const qglviewer::Camera* const camera)
{
    if (!applyManipulations) {
        setGrabsMouse(false);
        return;
    }
    Vector3 mousePos(x, y);
    Vector3 orig, dir;
    this->interactor_.computeMouseRay(mousePos, camera, orig, dir);

    this->interactor_.updateHover(state(), mousePos, camera);
    // state().visual.hovered = this->interactor_.pickHandle(state(), mousePos, camera);
    bool grabbed = state().visual.hovered != GizmoHandleType::None;
    this->setGrabsMouse(grabbed);
}

void Gizmo3D::mousePressEvent(QMouseEvent* const event, qglviewer::Camera* const camera)
{
    Vector3 mousePos(event->x(), event->y());
    this->interactor_.beginDrag(state(), mousePos, camera);
    qglviewer::MouseGrabber::mousePressEvent(event, camera);
}

void Gizmo3D::mouseMoveEvent(QMouseEvent* const event, qglviewer::Camera* const camera)
{
    if (isManipulated()) {
        Vector3 mousePos(event->x(), event->y());
        this->interactor_.updateDrag(state(), mousePos, camera);
        Q_EMIT this->pointModified();
    }
    qglviewer::MouseGrabber::mouseMoveEvent(event, camera);
}

void Gizmo3D::mouseReleaseEvent(QMouseEvent* const event, qglviewer::Camera* const camera)
{
    if (this->isManipulated()) {
        this->interactor_.endDrag(state());
        Q_EMIT this->pointReleased();
        Q_EMIT this->pointModified();
    }
    qglviewer::MouseGrabber::mouseReleaseEvent(event, camera);
}

void Gizmo3D::mouseDoubleClickEvent(QMouseEvent * const event, qglviewer::Camera * const camera)
{
    std::cout << "Mouse double-click" << std::endl;
    qglviewer::MouseGrabber::mouseDoubleClickEvent(event, camera);
}

void Gizmo3D::wheelEvent(QWheelEvent * const event, qglviewer::Camera * const camera)
{
    std::cout << "Mouse wheel event" << std::endl;
    qglviewer::MouseGrabber::wheelEvent(event, camera);
}

bool Gizmo3D::isManipulated() const
{
    return state().drag.active;
}
