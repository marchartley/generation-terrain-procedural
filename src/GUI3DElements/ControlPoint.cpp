#include "ControlPoint.h"


// #include "Utils/Collisions.h"


std::map<ControlPoint::GrabberState, Vector3> ControlPoint::default_GrabberStateColor = {
    {ControlPoint::GrabberState::HIDDEN, Vector3(.0f, .0f, .0f)},
    // {ControlPoint::GrabberState::INACTIVE, Vector3(.8f, .0f, .0f)},
    // {ControlPoint::GrabberState::ACTIVE, Vector3(.8f, .0f, .0f)},
    {ControlPoint::GrabberState::POSITIVE, Vector3(.2f, 1.f, .1f)},
    {ControlPoint::GrabberState::NEGATIVE, Vector3(1.f, .2f, .1f)},
    {ControlPoint::GrabberState::NEUTRAL, Vector3(.8f, .8f, .8f)},
    {ControlPoint::GrabberState::DEFAULT, Vector3(1.f, 0.f, 0.f)}
};

ControlPoint::ControlPoint()
    : ControlPoint(Vector3::origin)
{

}

ControlPoint::ControlPoint(const Vector3 &pos, float radius, GrabberState state, bool applyManipulations)
    : Gizmo3D(applyManipulations)
{
    this->stateColorMap = ControlPoint::default_GrabberStateColor;
    this->move(pos);
    this->setRadius(radius);
    this->setState(state);
}

Vector3 ControlPoint::getPosition() const {
    return this->state().transform.position;
}

Vector3 ControlPoint::getRotation() const
{
    return this->state().transform.rotation.toVector3();
}

void ControlPoint::allowAllAxisRotations(bool allow)
{
    state().capabilities.allowRotateX = allow;
    state().capabilities.allowRotateY = allow;
    state().capabilities.allowRotateZ = allow;
}

void ControlPoint::allowAllAxisTranslation(bool allow)
{
    state().capabilities.allowTranslateX = allow;
    state().capabilities.allowTranslateY = allow;
    state().capabilities.allowTranslateZ = allow;
}

void ControlPoint::display()
{
    if (state().transform.position.isValid())
        this->render();
}

void ControlPoint::hide()
{
    state().settings.visible = false;
    this->removeFromMouseGrabberPool();
}

void ControlPoint::show()
{
    state().settings.visible = true;
    this->addInMouseGrabberPool();
}

void ControlPoint::setVisible(bool visibility)
{
    if (visibility) show();
    else hide();
}

void ControlPoint::move(const Vector3 &newPos)
{
    state().transform.position = newPos;
    emitPointModified();
    // Q_EMIT this->pointModified();
}

void ControlPoint::setGrabberStateColor(GrabberState state, Vector3 color)
{
    this->setGrabberStateColor(std::map<GrabberState, Vector3>{{state, color}});
}

void ControlPoint::setGrabberStateColor(std::map<GrabberState, Vector3 > stateColorMap)
{
    for (auto [s, color] : stateColorMap) {
        this->stateColorMap[s] = color;
    }
}

void ControlPoint::setState(GrabberState newState)
{
    // std::cout << "ControlPoint3D::setState not implemented" << std::endl;
    state().visual.body_color = this->stateColorMap[newState];
}

void ControlPoint::setRadius(float newRadius)
{
    state().settings.radius = newRadius;
    state().settings.ringRadius = 1.5f * newRadius;
    state().settings.axisLength = 2.f * newRadius;
}

void ControlPoint::setConstraint(qglviewer::Constraint *constraint)
{
    std::cout << "ControlPoint3D::setConstraint not implemented" << std::endl;
}

// bool ControlPoint::grabsMouse() const
// {
    // std::cout << "ControlPoint::grabsMouse not implemented" << std::endl;
    // return false;
// }


void ControlPoint::setDisplayOnTop(bool enable)
{
    state().settings.displayOnTop = enable;
}

Vector3 ControlPoint::getFluidTranslation()
{
    std::cout << "ControlPoint3D::getFluidTranslation not implemented" << std::endl;
    return Vector3::invalid;
}





/*
std::shared_ptr<Shader> ControlPoint::base_shader = nullptr;
std::map<ControlPoint::GrabberState, std::vector<float>> ControlPoint::default_GrabberStateColor = {
    {ControlPoint::GrabberState::HIDDEN, {.0f, .0f, .0f, 0.f}},
    {ControlPoint::GrabberState::INACTIVE, {.3f, .0f, .0f, .5f}},
    {ControlPoint::GrabberState::ACTIVE, {.8f, .0f, .0f, .8f}},
    {ControlPoint::GrabberState::POSITIVE, {.2f, 1.f, .1f, .8f}},
    {ControlPoint::GrabberState::NEGATIVE, {1.f, .2f, .1f, .8f}},
    {ControlPoint::GrabberState::NEUTRAL, {.8f, .8f, .8f, .8f}},
    };

ControlPoint::ControlPoint()
    : ControlPoint(Vector3())
{
    this->mesh.hide();
    this->rotationMeshes.show();
    this->translationMeshes.show();
    this->allowAllAxisRotations(false);
    this->allowAllAxisTranslation(false);
}

ControlPoint::ControlPoint(const Vector3& pos, float radius, GrabberState state, bool useTheManipulatedFrame)
    : state(state), useManipFrame(useTheManipulatedFrame), radius(radius)
{
    this->mesh = Mesh((ControlPoint::base_shader ? std::make_shared<Shader>(ControlPoint::base_shader->vertexShaderFilename, ControlPoint::base_shader->fragmentShaderFilename, ControlPoint::base_shader->geometryShaderFilename) : nullptr), true);
    this->move(pos);
    this->prevPosition = pos;
    this->GrabberStateColor = ControlPoint::default_GrabberStateColor;
    this->currentlyManipulated = false;
    if (!useManipFrame) {
        this->removeFromMouseGrabberPool();
    } else if (!this->isInMouseGrabberPool()) {
        this->addInMouseGrabberPool();
    }
    this->allowAllAxisRotations(false);
    this->allowAllAxisTranslation(false);
    this->rotationMeshes.show();
    this->translationMeshes.show();

    QObject::connect(this, &qglviewer::ManipulatedFrame::modified, this, [=](){
        Q_EMIT ControlPoint::pointModified();
        if ((this->prevPosition - this->getPosition()).norm2() > 1.0) {
            this->prevPosition = this->getPosition();
//        if (this->positionsHistory.empty() || this->positionsHistory.back() != this->prevPosition) {
            this->positionsHistory.push_back(prevPosition);
            if (this->positionsHistory.size() > 10) {
                this->positionsHistory.erase(positionsHistory.begin(), std::max(positionsHistory.end() - 10, positionsHistory.begin()));
            }
        }
        this->updateStateDependingOnManipFrame();
    });
}

ControlPoint::~ControlPoint()
{
    this->removeFromMouseGrabberPool();
}

void ControlPoint::setState(GrabberState newState)
{
    this->state = newState;
}

void ControlPoint::setVisible(bool visibility)
{
    this->visible = visibility;
}

void ControlPoint::updateStateDependingOnManipFrame()
{
    if (this->useManipFrame && this->visible)
        this->setState(this->isManipulated() ? ACTIVE : INACTIVE);
}

void ControlPoint::checkIfGrabsMouse(int x, int y, const qglviewer::Camera * const cam)
{
    if (this->isManipulated())
    {
        setGrabsMouse(true);
        return;
    }
    if (!this->visible || this->mesh.isHidden()) {
        setGrabsMouse(false);
        return;
    }

    qglviewer::Vec orig, dir;
    cam->convertClickToLine(QPoint(x,y), orig, dir);
    Vector3 rayOrigin = Vector3(orig), rayDir = Vector3(dir);

    auto [rotationAxis, intersectionPoint] = this->mouseOnRotationCircle(rayOrigin, rayDir);

    if (this->mouseOnCentralSphere(rayOrigin, rayDir)) {
        setGrabsMouse(true);
        this->interactionState.currentAxis = ControlPointAxis::NONE;
        this->interactionState.isApplyingFreeMove = true;
        this->interactionState.isApplyingRotation = false;
        this->interactionState.isApplyingTranslation = false;
    }
    else if (this->mouseOnTranslationArrow(rayOrigin, rayDir)) {
        setGrabsMouse(true);
        this->interactionState.isApplyingFreeMove = false;
        this->interactionState.isApplyingRotation = false;
        this->interactionState.isApplyingTranslation = true;
    } else if (rotationAxis != ControlPointAxis::NONE) {
        setGrabsMouse(true);
        this->interactionState.isApplyingFreeMove = false;
        this->interactionState.isApplyingRotation = true;
        this->interactionState.isApplyingTranslation = false;
        this->interactionState.currentMousePosOnAction = intersectionPoint;
        this->interactionState.currentAxis = rotationAxis;
    } else {
        setGrabsMouse(false);
        this->interactionState.isApplyingFreeMove = false;
        this->interactionState.isApplyingRotation = false;
        this->interactionState.isApplyingTranslation = false;
    }


    // Constraints :
    if (this->default_constraint == nullptr)
        this->default_constraint = new qglviewer::WorldConstraint();
    auto& constraint = this->default_constraint;

    if (this->interactionState.isApplyingFreeMove) {
        constraint->setTranslationConstraintType(qglviewer::AxisPlaneConstraint::FREE);
        constraint->setRotationConstraintType(qglviewer::AxisPlaneConstraint::FREE);
    } else if (this->interactionState.isApplyingRotation) {
        constraint->setRotationConstraintType(qglviewer::AxisPlaneConstraint::AXIS);
        if (this->interactionState.currentAxis == ControlPointAxis::X)
            constraint->setRotationConstraintDirection(Vector3(1.0, 0.0, 0.0));
        else if (this->interactionState.currentAxis == ControlPointAxis::Y)
            constraint->setRotationConstraintDirection(Vector3(0.0, 1.0, 0.0));
        else if (this->interactionState.currentAxis == ControlPointAxis::Z)
            constraint->setRotationConstraintDirection(Vector3(0.0, 0.0, 1.0));
    } else if (this->interactionState.isApplyingTranslation) {
        constraint->setTranslationConstraintType(qglviewer::AxisPlaneConstraint::AXIS);
        if (this->interactionState.currentAxis == ControlPointAxis::X)
            constraint->setTranslationConstraintDirection(Vector3(1.0, 0.0, 0.0));
        else if (this->interactionState.currentAxis == ControlPointAxis::Y)
            constraint->setTranslationConstraintDirection(Vector3(0.0, 1.0, 0.0));
        else if (this->interactionState.currentAxis == ControlPointAxis::Z)
            constraint->setTranslationConstraintDirection(Vector3(0.0, 0.0, 1.0));
    }
    if (this->custom_constraint == nullptr) {
        this->setConstraint(constraint);
    } else {
        this->setConstraint(this->custom_constraint);
    }
}

void ControlPoint::onUpdate(std::function<void ()> func)
{
    this->onUpdateCallback = func;
}

void ControlPoint::afterUpdate(std::function<void ()> func)
{
    this->afterUpdateCallback = func;
}

void ControlPoint::updateInteractionState()
{

}

void ControlPoint::updateGeometryIfNeeded()
{
    if (!geometryDirty) return;

    Sphere s(this->radius, getPosition(), 10, 10);
    s.buildVerticesFlat();
    mesh.fromArray(s.mesh.vertexArrayFloat);
    mesh.update();

    geometryDirty = false;
}

void ControlPoint::updateSphere()
{
    if (this->useManipFrame && (this->isManipulated() == true && this->currentlyManipulated == false)) {
        this->initialPosition = this->getPosition();
//        this->initialRotation = this->getRotation();
    }

    if(this->onUpdateCallback) { // && this->useManipFrame && this->isManipulated())
        if (prevPosition != getPosition()) {
            this->prevPosition = getPosition();
            this->onUpdateCallback();
        }
    }
    if (this->afterUpdateCallback) {
        if (this->useManipFrame && (this->isManipulated() == false && this->currentlyManipulated == true)) {
            this->afterUpdateCallback();
        }
    }
    if (this->useManipFrame && (this->isManipulated() == false && this->currentlyManipulated == true)) {
        Q_EMIT this->pointReleased();
        Vector3 translation = this->getPosition() - initialPosition;
        Vector3 rotation = this->getRotation();
        if (translation.norm2() > 0) {
            Q_EMIT this->translationApplied(translation);
        }
        if (rotation.abs().maxComp() > 0) {
            Q_EMIT this->rotationApplied(rotation); // - initialRotation);
            this->setRotation(qglviewer::Quaternion()); // Back to identity (?)
        }
    }
    this->currentlyManipulated = this->isManipulated();
//    this->manipFrame.setPosition(this->getPosition());

    Sphere s(this->radius, this->getPosition(), 6, 6);
    s.buildVerticesFlat();
    this->mesh.fromArray(s.mesh.vertexArrayFloat);
    this->mesh.update();

    this->arrowSize = 3 * this->radius;
    this->circleRadius = 2 * this->radius;
}

void ControlPoint::display()
{
    if (!this->visible) return;

    this->updateGeometryIfNeeded();

    GLboolean m_origin_blend, m_origin_depth, m_origin_cull;
    glGetBooleanv(GL_BLEND, &m_origin_blend);
    glGetBooleanv(GL_DEPTH_TEST,&m_origin_depth);
    glGetBooleanv(GL_CULL_FACE, &m_origin_cull);
    if (this->displayOnTop) {
        glEnable(GL_BLEND);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);
    }
    float controlAxisSizeUnselected = 2.f;
    float controlAxisSizeSelected = 4.f;
    if (this->useManipFrame)
        this->setState(this->isManipulated() ? ACTIVE : INACTIVE);

    this->mesh.shader->setVector("color", ControlPoint::GrabberStateColor[this->state]);

    if (this->translationMeshes.shader != nullptr) {
        // Display X (red)
        if (this->allowedTranslations[ControlPointAxis::X]) {
            this->translationMeshes.shader->setVector("color", std::vector<float>({1.0, 0.0, 0.0, 1.0}));
            this->translationMeshes.fromArray({this->getPosition() - Vector3(1.0, 0.0, 0.0) * arrowSize, this->getPosition() + Vector3(1.0, 0.0, 0.0) * arrowSize});
            this->translationMeshes.display(GL_LINES, (this->interactionState.isApplyingTranslation && this->interactionState.currentAxis == ControlPointAxis::X ? controlAxisSizeSelected : controlAxisSizeUnselected));
        }
        // Display Y (green)
        if (this->allowedTranslations[ControlPointAxis::Y]) {
            this->translationMeshes.shader->setVector("color", std::vector<float>({0.0, 1.0, 0.0, 1.0}));
            this->translationMeshes.fromArray({this->getPosition() - Vector3(0.0, 1.0, 0.0) * arrowSize, this->getPosition() + Vector3(0.0, 1.0, 0.0) * arrowSize});
            this->translationMeshes.display(GL_LINES, (this->interactionState.isApplyingTranslation && this->interactionState.currentAxis == ControlPointAxis::Y ? controlAxisSizeSelected : controlAxisSizeUnselected));
        }
        // Display Z (blue)
        if (this->allowedTranslations[ControlPointAxis::Z]) {
            this->translationMeshes.shader->setVector("color", std::vector<float>({0.0, 0.0, 1.0, 1.0}));
            this->translationMeshes.fromArray({this->getPosition() - Vector3(0.0, 0.0, 1.0) * arrowSize, this->getPosition() + Vector3(0.0, 0.0, 1.0) * arrowSize});
            this->translationMeshes.display(GL_LINES, (this->interactionState.isApplyingTranslation && this->interactionState.currentAxis == ControlPointAxis::Z ? controlAxisSizeSelected : controlAxisSizeUnselected));
        }
    } else if (this->mesh.shader != nullptr ){
        this->translationMeshes.shader = std::make_shared<Shader>(*this->mesh.shader);
    }
    if (this->rotationMeshes.shader != nullptr) {
        // Display X (red)
        if (this->allowedRotations[ControlPointAxis::X]) {
            this->rotationMeshes.shader->setVector("color", std::vector<float>({1.0, 0.0, 0.0, 1.0}));
            this->rotationMeshes.fromArray(computeCircle(ControlPointAxis::X));
            this->rotationMeshes.display(GL_LINES, (this->interactionState.isApplyingRotation && this->interactionState.currentAxis == ControlPointAxis::X ? controlAxisSizeSelected : controlAxisSizeUnselected));
            if (this->interactionState.isApplyingRotation && this->interactionState.currentAxis == ControlPointAxis::X) {
//                    this->rotationHelperSphere.translate(Vector3(1, 0, 0));
                this->rotationHelperSphere.display();
            }
        }
        // Display Y (green)
        if (this->allowedRotations[ControlPointAxis::Y]) {
            this->rotationMeshes.shader->setVector("color", std::vector<float>({0.0, 1.0, 0.0, 1.0}));
            this->rotationMeshes.fromArray(computeCircle(ControlPointAxis::Y));
            this->rotationMeshes.display(GL_LINES, (this->interactionState.isApplyingRotation && this->interactionState.currentAxis == ControlPointAxis::Y ? controlAxisSizeSelected : controlAxisSizeUnselected));
            if (this->interactionState.isApplyingRotation && this->interactionState.currentAxis == ControlPointAxis::Y) {
//                    this->rotationHelperSphere.translate(Vector3(1, 0, 0));
                this->rotationHelperSphere.display();
            }
        }
        // Display Z (blue)
        if (this->allowedRotations[ControlPointAxis::Z]) {
            this->rotationMeshes.shader->setVector("color", std::vector<float>({0.0, 0.0, 1.0, 1.0}));
            this->rotationMeshes.fromArray(computeCircle(ControlPointAxis::Z));
            this->rotationMeshes.display(GL_LINES, (this->interactionState.isApplyingRotation && this->interactionState.currentAxis == ControlPointAxis::Z ? controlAxisSizeSelected : controlAxisSizeUnselected));
            if (this->interactionState.isApplyingRotation && this->interactionState.currentAxis == ControlPointAxis::Z) {
//                    this->rotationHelperSphere.translate(Vector3(1, 0, 0));
                this->rotationHelperSphere.display();
            }
        }
    } else if (this->mesh.shader != nullptr ){
        this->rotationMeshes.shader = std::make_shared<Shader>(*this->mesh.shader);
        this->rotationHelperSphere.shareShader(rotationMeshes.shader);
    }
    this->updateSphere();
    this->mesh.display();

    if (m_origin_blend == GL_TRUE) glEnable(GL_BLEND);
    else glDisable(GL_BLEND);
    if (m_origin_depth == GL_TRUE) glEnable(GL_DEPTH_TEST);
    else glDisable(GL_DEPTH_TEST);
    if (m_origin_cull == GL_TRUE) glEnable(GL_CULL_FACE);
    else glDisable(GL_CULL_FACE);
}

void ControlPoint::hide()
{
    this->removeFromMouseGrabberPool();
    this->mesh.hide();
    this->translationMeshes.hide();
    this->rotationMeshes.hide();
    this->rotationHelperSphere.hide();
//    CustomInteractiveObject::hide();
}

void ControlPoint::show()
{
    if (useManipFrame)
        this->addInMouseGrabberPool();
    this->mesh.show();
    this->translationMeshes.show();
    this->rotationMeshes.show();
    this->rotationHelperSphere.show();
//    CustomInteractiveObject::show();
}

void ControlPoint::setGrabberStateColor(std::map<GrabberState, std::vector<float> > stateColorMap)
{
    for (auto& tuple : stateColorMap) {
        this->setGrabberStateColor(std::get<0>(tuple), std::get<1>(tuple));
    }
}

void ControlPoint::setGrabberStateColor(GrabberState state, std::vector<float> color)
{
    this->GrabberStateColor[state] = color;
}

void ControlPoint::setPosition(const Vector3& newPos)
{
    this->setPosition(newPos.x(), newPos.y(), newPos.z());
}

void ControlPoint::setPosition(float x, float y, float z)
{
    this->currentPosition = Vector3(x, y, z);
    ManipulatedFrame::setPosition(x, y, z);
}

Vector3 ControlPoint::getRotation() const
{
    if (!this->interactionState.pressedPosBeforeAction.isValid())
        return Vector3(0.f, 0.f, 0.f);
    Vector3 rotation = this->interactionState.pressedPosBeforeAction.getAllAnglesWith(this->interactionState.currentMousePosOnAction);
    return rotation;
}

Vector3 ControlPoint::getPosition() const
{
    return this->position();
    // return this->currentPosition;
}

Vector3 ControlPoint::getFluidTranslation() const {
    if (positionsHistory.empty()) return Vector3();
    return (this->getPosition() - this->positionsHistory.front()).normalize(); }

Vector3 ControlPoint::getLastMovement() const
{
    return (this->getPosition() - this->prevPosition).normalize();
}

void ControlPoint::mousePressEvent(QMouseEvent * const event, qglviewer::Camera * const cam)
{
    std::cout << "ControlPoint::mousePressEvent triggered" << std::endl;
    this->interactionState.pressedPosBeforeAction.setValid(false);
    if (this->grabsMouse()) {
        if (this->interactionState.isApplyingRotation) {
            this->interactionState.pressedPosBeforeAction = interactionState.currentMousePosOnAction;
            Sphere s(1.f, this->getPosition() + circleRadius * this->interactionState.pressedPosBeforeAction.normalized(), 6, 6);
            s.buildVerticesFlat();
            this->rotationHelperSphere.fromArray(s.mesh.vertexArray);
            this->startAction(QGLViewer::ROTATE);
        } else if (this->interactionState.isApplyingTranslation) {
            this->startAction(QGLViewer::TRANSLATE);
        } else if (this->interactionState.isApplyingFreeMove) {
            this->startAction(QGLViewer::TRANSLATE); // force translation
        }
    }
    qglviewer::ManipulatedFrame::mousePressEvent(event, cam);
}

void ControlPoint::mouseReleaseEvent(QMouseEvent * const event, qglviewer::Camera * const cam)
{
    this->currentlyManipulated = false;
    qglviewer::ManipulatedFrame::mouseReleaseEvent(event, cam);
}

void ControlPoint::mouseMoveEvent(QMouseEvent * const event, qglviewer::Camera * const cam)
{
    if (this->interactionState.isApplyingRotation) {
        Sphere s(1.f, this->getPosition() + circleRadius * this->interactionState.currentMousePosOnAction.normalized(), 6, 6);
        s.buildVerticesFlat();
        this->rotationHelperSphere.fromArray(s.mesh.vertexArray);
    }
    event->accept();
    qglviewer::ManipulatedFrame::mouseMoveEvent(event, cam);
}

void ControlPoint::wheelEvent(QWheelEvent * const event, qglviewer::Camera * const camera)
{
    setRadius(this->radius - event->angleDelta().y()/10.f);
    this->startAction(QGLViewer::MouseAction::NO_MOUSE_ACTION);
    qglviewer::ManipulatedFrame::wheelEvent(event, camera);
}

void ControlPoint::setRadius(float newRadius)
{
    if (minSphereRadius >= 0.0f) newRadius = std::max(minSphereRadius, newRadius);
    if (maxSphereRadius >= 0.0f) newRadius = std::min(maxSphereRadius, newRadius);

    if (radius != newRadius) {
        radius = newRadius;
        arrowSize = 3.0f * radius;
        circleRadius = 2.0f * radius;
        geometryDirty = true;
    }
}

void ControlPoint::allowAllAxisTranslation(bool allow)
{
    this->allowedTranslations[ControlPointAxis::X] = allow;
    this->allowedTranslations[ControlPointAxis::Y] = allow;
    this->allowedTranslations[ControlPointAxis::Z] = allow;
}

void ControlPoint::allowAllAxisRotations(bool allow)
{
    this->allowedRotations[ControlPointAxis::X] = allow;
    this->allowedRotations[ControlPointAxis::Y] = allow;
    this->allowedRotations[ControlPointAxis::Z] = allow;
}

void ControlPoint::setConstraint(qglviewer::Constraint * const constraint)
{
    this->custom_constraint = constraint;
    ManipulatedFrame::setConstraint(constraint);
}

std::vector<Vector3> ControlPoint::computeCircle(ControlPointAxis axis)
{
    std::vector<Vector3> points;
    for (int i = 0; i <= 360; i += 5) {
        float angle = i * PI / 180.f;
        float nextAngle = (i + 5) * PI / 180.f;
        if (axis == ControlPointAxis::X) {
            points.push_back(this->getPosition() + Vector3(0.0, 1.0, 0.0).rotate(angle, 0, 0) * this->circleRadius);
            points.push_back(this->getPosition() + Vector3(0.0, 1.0, 0.0).rotate(nextAngle, 0, 0) * this->circleRadius);
        } else if (axis == ControlPointAxis::Y) {
            points.push_back(this->getPosition() + Vector3(1.0, 0.0, 0.0).rotate(0, angle, 0) * this->circleRadius);
            points.push_back(this->getPosition() + Vector3(1.0, 0.0, 0.0).rotate(0, nextAngle, 0) * this->circleRadius);
        } else if (axis == ControlPointAxis::Z) {
            points.push_back(this->getPosition() + Vector3(1.0, 0.0, 0.0).rotate(0, 0, angle) * this->circleRadius);
            points.push_back(this->getPosition() + Vector3(1.0, 0.0, 0.0).rotate(0, 0, nextAngle) * this->circleRadius);
        }
    }
    return points;
}

Vector3 ControlPoint::getIntersectionWithPlane(const Vector3& rayOrigin, const Vector3& rayDir, ControlPointAxis axis)
{
    return Collision::intersectionRayPlane(rayOrigin, rayDir, this->getPosition(), Vector3(
                                               (axis == ControlPointAxis::X ? 1.f :  0.f),
                                               (axis == ControlPointAxis::Y ? 1.f :  0.f),
                                               (axis == ControlPointAxis::Z ? 1.f :  0.f)));
}

bool ControlPoint::mouseOnCentralSphere(const Vector3& rayOrigin, const Vector3& rayDir)
{
    return Collision::intersectionRaySphere(rayOrigin, rayDir, this->getPosition(), this->radius).isValid();
}

bool ControlPoint::mouseOnTranslationArrow(const Vector3& rayOrigin, const Vector3& rayDir)
{
    float distToCam = (this->getPosition() - rayOrigin).norm();
    float tolerance = this->radius/(650.f / distToCam); // Magic number, I know...

    // X-axis
    if (this->allowedTranslations[ControlPointAxis::X]) {
        float distance = Collision::shortestDistanceBetweenSegments(rayOrigin, rayOrigin + rayDir * 1000.f,
                                                                    this->getPosition() - Vector3(1.0, 0.0, 0.0) * arrowSize,
                                                                    this->getPosition() + Vector3(1.0, 0.0, 0.0) * arrowSize);
        if(distance < tolerance) {
            this->interactionState.currentAxis = ControlPointAxis::X;
            return true;
        }
    }
    // Y-axis
    if (this->allowedTranslations[ControlPointAxis::Y]) {
        float distance = Collision::shortestDistanceBetweenSegments(rayOrigin, rayOrigin + rayDir * 1000.f,
                                        this->getPosition() - Vector3(0.0, 1.0, 0.0) * arrowSize,
                                        this->getPosition() + Vector3(0.0, 1.0, 0.0) * arrowSize);
        if (distance < tolerance) {
            this->interactionState.currentAxis = ControlPointAxis::Y;
            return true;
        }
    }
    // Z-axis
    if (this->allowedTranslations[ControlPointAxis::Z]) {
        float distance = Collision::shortestDistanceBetweenSegments(rayOrigin, rayOrigin + rayDir * 1000.f,
                                        this->getPosition() - Vector3(0.0, 0.0, 1.0) * arrowSize,
                                        this->getPosition() + Vector3(0.0, 0.0, 1.0) * arrowSize);
        if (distance < tolerance) {
            this->interactionState.currentAxis = ControlPointAxis::Z;
            return true;
        }
    }
    return false;
}

std::pair<ControlPointAxis, Vector3> ControlPoint::mouseOnRotationCircle(const Vector3& rayOrigin, const Vector3& rayDir)
{
//    float circleRadiusSq = this->circleRadius * this->circleRadius;
    float tolerence = this->radius / 5.f;
    float minCircleRadiusSq = (this->circleRadius - tolerence) * (this->circleRadius - tolerence);
    float maxCircleRadiusSq = (this->circleRadius + tolerence) * (this->circleRadius + tolerence);
    Vector3 intersection;
    Vector3 bestIntersection;
    float distanceToCamSq = std::numeric_limits<float>::max();
    float distanceToCenterSq;
    ControlPointAxis bestAxis = ControlPointAxis::NONE;

    // X-axis
    if (this->allowedRotations[ControlPointAxis::X]) {
        intersection = this->getIntersectionWithPlane(rayOrigin, rayDir, ControlPointAxis::X);
        if (intersection.isValid()) {
            distanceToCenterSq = (intersection - this->getPosition()).norm2();
            if (minCircleRadiusSq < distanceToCenterSq && distanceToCenterSq < maxCircleRadiusSq) {
                if ((intersection - rayOrigin).norm2() < distanceToCamSq) {
                    distanceToCamSq = (intersection - rayOrigin).norm2();
                    bestAxis = ControlPointAxis::X;
                    bestIntersection = intersection;
                }
            }
        }
    }
    // Y-axis
    if (this->allowedRotations[ControlPointAxis::Y]) {
        intersection = this->getIntersectionWithPlane(rayOrigin, rayDir, ControlPointAxis::Y);
        if (intersection.isValid()) {
            distanceToCenterSq = (intersection - this->getPosition()).norm2();
            if (minCircleRadiusSq < distanceToCenterSq && distanceToCenterSq < maxCircleRadiusSq) {
                if ((intersection - rayOrigin).norm2() < distanceToCamSq) {
                    distanceToCamSq = (intersection - rayOrigin).norm2();
                    bestAxis = ControlPointAxis::Y;
                    bestIntersection = intersection;
                }
            }
        }
    }
    // Z-axis
    if (this->allowedRotations[ControlPointAxis::Z]) {
        intersection = this->getIntersectionWithPlane(rayOrigin, rayDir, ControlPointAxis::Z);
        if (intersection.isValid()) {
            distanceToCenterSq = (intersection - this->getPosition()).norm2();
            if (minCircleRadiusSq < distanceToCenterSq && distanceToCenterSq < maxCircleRadiusSq) {
                if ((intersection - rayOrigin).norm2() < distanceToCamSq) {
                    distanceToCamSq = (intersection - rayOrigin).norm2();
                    bestAxis = ControlPointAxis::Z;
                    bestIntersection = intersection;
                }
            }
        }
    }
    return {bestAxis, bestIntersection - this->getPosition()};
}

void ControlPoint::move(const Vector3& newPos)
{
    if (!this->isManipulated()) {
        this->setPosition(newPos.x(), newPos.y(), newPos.z());
        this->updateSphere();
    }
    this->currentPosition = newPos;
}
*/



CustomConstraint::CustomConstraint()
{
    this->constraint = new qglviewer::WorldConstraint();
}

void CustomConstraint::constrainTranslation(qglviewer::Vec& t, qglviewer::Frame* const fr) {
    this->constraint->constrainTranslation(t, fr);
}

void CustomConstraint::constrainRotation(qglviewer::Quaternion &q, qglviewer::Frame * const fr)
{
    this->constraint->constrainRotation(q, fr);
}
