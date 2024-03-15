#include "ControlPoint.h"

#include "Utils/Collisions.h"


std::shared_ptr<Shader> ControlPoint::base_shader = nullptr;
std::map<GrabberState, std::vector<float>> ControlPoint::default_GrabberStateColor = {
    {GrabberState::HIDDEN, {.0f, .0f, .0f, 0.f}},
    {GrabberState::INACTIVE, {.3f, .0f, .0f, .5f}},
    {GrabberState::ACTIVE, {.8f, .0f, .0f, .8f}},
    {GrabberState::POSITIVE, {.2f, 1.f, .1f, .8f}},
    {GrabberState::NEGATIVE, {1.f, .2f, .1f, 8.f}},
    {GrabberState::NEUTRAL, {.8f, .8f, .8f, .8f}},
};

/*
ControlPoint::ControlPoint()
    : ControlPoint(Vector3())
{}

ControlPoint::ControlPoint(const Vector3& pos, float radius, GrabberState state, bool useTheManipulatedFrame)
    : position(pos), radius(radius), currentState(state)
{
    this->allowAllAxisTranslation(true);
    this->allowAllAxisRotations(true);
    this->allowAllAxisScaling(true);
//    translationMesh.show();
//    rotationMesh.show();
//    scalingMesh.show();
    this->addInMouseGrabberPool();

    // Need to take into account "custom_constraint" for translations on specific axis
}

ControlPoint::~ControlPoint()
{
    this->removeFromMouseGrabberPool();
}

void ControlPoint::display()
{
//    translationMesh.display();
//    rotationMesh.display();
//    scalingMesh.display();
}

void ControlPoint::move(const Vector3& newPos)
{
    this->setPosition(newPos);
}

Vector3 ControlPoint::getPosition() const
{
    return this->position;
}

void ControlPoint::allowAllAxisTranslation(bool allowed)
{
    this->allowedTranslations[X] = allowed;
    this->allowedTranslations[Y] = allowed;
    this->allowedTranslations[Z] = allowed;
}

void ControlPoint::allowAllAxisRotations(bool allowed)
{
    this->allowedRotations[X] = allowed;
    this->allowedRotations[Y] = allowed;
    this->allowedRotations[Z] = allowed;
}

void ControlPoint::allowAllAxisScaling(bool allowed)
{
    this->allowedScaling[X] = allowed;
    this->allowedScaling[Y] = allowed;
    this->allowedScaling[Z] = allowed;
}

void ControlPoint::mousePressEvent(QMouseEvent * const event, qglviewer::Camera * const cam)
{
    if (this->grabsMouse()) {
        std::cout << "Grabbed!" << std::endl;
    }
    // Check if visible
    // Check grabMouse
    // Get action
    // Get axis
    // Save "pressedPos"
}

void ControlPoint::mouseReleaseEvent(QMouseEvent * const event, qglviewer::Camera * const cam)
{
//    Q_EMIT this->translationApplied();
//    Q_EMIT this->rotationApplied();
    Q_EMIT this->released();
}

void ControlPoint::mouseMoveEvent(QMouseEvent * const event, qglviewer::Camera * const cam)
{
    Q_EMIT this->modified();
}

void ControlPoint::wheelEvent(QWheelEvent * const event, qglviewer::Camera * const camera)
{
    this->setRadius(radius - event->angleDelta().y()/10.f);
    // Update meshes
}

void ControlPoint::hide()
{
//    translationMesh.hide();
//    rotationMesh.hide();
//    scalingMesh.hide();
}

void ControlPoint::show()
{
//    translationMesh.show();
//    rotationMesh.show();
//    scalingMesh.show();
}

void ControlPoint::setPosition(const Vector3& newPosition)
{
    this->position = newPosition;
    // Update meshes?
}

void ControlPoint::setRadius(float newRadius)
{
    this->radius = std::max(newRadius, .5f);
}

void ControlPoint::setState(GrabberState newState)
{
    // Update meshes
}

Vector3 ControlPoint::getFluidTranslation() const
{
    return this->getCurrentTranslation();
}

Vector3 ControlPoint::getFluidRotation() const
{
    return this->getCurrentRotation();
}

Vector3 ControlPoint::getCurrentTranslation() const
{

}

Vector3 ControlPoint::getCurrentRotation() const
{

}

void ControlPoint::checkIfGrabsMouse(int x, int y, const qglviewer::Camera * const cam)
{
    // Check intersection with all widgets
}

Vector3 ControlPoint::intersectionWithTranslationWidget(const Vector3& rayOrigin, const Vector3& rayDir)
{
    // Check mouse-segments with X, Y, Z. return closest to cam
    Vector3 intersectX = this->getIntersectionWithTranslationAxis(rayOrigin, rayDir, X);
    Vector3 intersectY = this->getIntersectionWithTranslationAxis(rayOrigin, rayDir, Y);
    Vector3 intersectZ = this->getIntersectionWithTranslationAxis(rayOrigin, rayDir, Z);
    if (!intersectX.isValid() && !intersectY.isValid() && !intersectZ.isValid())
        return Vector3(false);

    float distX = (intersectX.isValid() ? (intersectX - rayOrigin).norm2() : std::numeric_limits<float>::max());
    float distY = (intersectY.isValid() ? (intersectY - rayOrigin).norm2() : std::numeric_limits<float>::max());
    float distZ = (intersectZ.isValid() ? (intersectZ - rayOrigin).norm2() : std::numeric_limits<float>::max());

    float minDist = std::min({distX, distY, distZ});
    if (distX == minDist) return intersectX;
    else if (distY == minDist) return intersectY;
    else return intersectZ;
}

Vector3 ControlPoint::intersectionWithRotationWidget(const Vector3& rayOrigin, const Vector3& rayDir)
{
    // Check mouse-plane
    Vector3 intersectX = this->getIntersectionWithRotationAxis(rayOrigin, rayDir, X);
    Vector3 intersectY = this->getIntersectionWithRotationAxis(rayOrigin, rayDir, Y);
    Vector3 intersectZ = this->getIntersectionWithRotationAxis(rayOrigin, rayDir, Z);
    if (!intersectX.isValid() && !intersectY.isValid() && !intersectZ.isValid())
        return Vector3(false);

    float distX = (intersectX.isValid() ? (intersectX - rayOrigin).norm2() : std::numeric_limits<float>::max());
    float distY = (intersectY.isValid() ? (intersectY - rayOrigin).norm2() : std::numeric_limits<float>::max());
    float distZ = (intersectZ.isValid() ? (intersectZ - rayOrigin).norm2() : std::numeric_limits<float>::max());

    float minDist = std::min({distX, distY, distZ});
    if (distX == minDist) return intersectX;
    else if (distY == minDist) return intersectY;
    else return intersectZ;
}

Vector3 ControlPoint::intersectionWithScalingWidget(const Vector3& rayOrigin, const Vector3& rayDir)
{
    // Check mouse-box
    Vector3 intersectX = this->getIntersectionWithScalingAxis(rayOrigin, rayDir, X);
    Vector3 intersectY = this->getIntersectionWithScalingAxis(rayOrigin, rayDir, Y);
    Vector3 intersectZ = this->getIntersectionWithScalingAxis(rayOrigin, rayDir, Z);
    if (!intersectX.isValid() && !intersectY.isValid() && !intersectZ.isValid())
        return Vector3(false);

    float distX = (intersectX.isValid() ? (intersectX - rayOrigin).norm2() : std::numeric_limits<float>::max());
    float distY = (intersectY.isValid() ? (intersectY - rayOrigin).norm2() : std::numeric_limits<float>::max());
    float distZ = (intersectZ.isValid() ? (intersectZ - rayOrigin).norm2() : std::numeric_limits<float>::max());

    float minDist = std::min({distX, distY, distZ});
    if (distX == minDist) return intersectX;
    else if (distY == minDist) return intersectY;
    else return intersectZ;
}

ControlPoint::ControlPointAction ControlPoint::hoveredAction(const Vector3& rayOrigin, const Vector3& rayDir)
{
    // Get intersection position with all widgets
    // Select closest to cam
}

ControlPoint::Axis ControlPoint::hoveredAxis(const Vector3& rayOrigin, const Vector3& rayDir)
{
    // From hovered action, check each component individually
}

Vector3 ControlPoint::getIntersectionWithTranslationAxis(const Vector3& rayOrigin, const Vector3& rayDir, Axis axis)
{
    float tolerence = this->radius * .2f;
    float arrowSize = this->radius * 2.f;
    Vector3 widgetAxis = Vector3((axis == X ? 1.f : 0.f), (axis == Y ? 1.f : 0.f), (axis == Z ? 1.f : 0.f)) * arrowSize;

    Vector3 intersection = Collision::intersectionBetweenTwoSegments(
                rayOrigin,
                rayOrigin + rayDir * std::max(arrowSize * 2.f, (this->getPosition() - rayOrigin).norm() * 2.f),
                this->getPosition() - widgetAxis,
                this->getPosition() + widgetAxis
                );
    if (intersection.isValid()) {
        Vector3 projection = Collision::projectPointOnSegment(intersection, this->getPosition() - widgetAxis, this->getPosition() + widgetAxis);
        if ((projection - intersection).norm2() > tolerence * tolerence)
            return Vector3(false);
        else
            return projection;
    } else {
        return Vector3(false);
    }
}

Vector3 ControlPoint::getIntersectionWithRotationAxis(const Vector3& rayOrigin, const Vector3& rayDir, Axis axis)
{
    float tolerence = this->radius * .2f;
    float circleSize = this->radius * 2.f;
    Vector3 widgetAxis = Vector3((axis == X ? 1.f : 0.f), (axis == Y ? 1.f : 0.f), (axis == Z ? 1.f : 0.f));

    Vector3 intersection = Collision::intersectionRayPlane(
                rayOrigin,
                rayOrigin + rayDir * std::max(circleSize * 2.f, (this->getPosition() - rayOrigin).norm() * 2.f),
                this->getPosition(),
                widgetAxis
                );
    if (intersection.isValid()) {
        Vector3 projection = Collision::projectPointOnSphere(intersection, this->getPosition(), circleSize);
        if ((projection - intersection).norm2() > tolerence * tolerence)
            return Vector3(false);
        else
            return projection;

    } else {
        return Vector3(false);
    }
}

Vector3 ControlPoint::getIntersectionWithScalingAxis(const Vector3& rayOrigin, const Vector3& rayDir, Axis axis)
{
    float tolerence = this->radius * .2f;
    float boxSize = this->radius * 2.f;
    float boxDistance = this->radius * 2.f;
    Vector3 widgetAxis = Vector3((axis == X ? 1.f : 0.f), (axis == Y ? 1.f : 0.f), (axis == Z ? 1.f : 0.f)) * boxDistance;

    Vector3 intersection = Collision::intersectionRayAABBox(
                rayOrigin,
                rayOrigin + rayDir * std::max(boxSize * 2.f, (this->getPosition() - rayOrigin).norm() * 2.f),
                this->getPosition() + widgetAxis - Vector3(boxSize, boxSize, boxSize) * .5f,
                this->getPosition() + widgetAxis + Vector3(boxSize, boxSize, boxSize) * .5f
                );
    if (intersection.isValid()) {
        return intersection;

    } else {
        return Vector3(false);
    }
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


*/


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
    : state(state), useManipFrame(useTheManipulatedFrame), shape(radius, getPosition(), 10, 10), radius(radius)
{
    this->mesh = Mesh((ControlPoint::base_shader ? std::make_shared<Shader>(*ControlPoint::base_shader) : nullptr), true);
    this->move(pos);
    this->stillOnInitialState = false; // true;
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

void ControlPoint::updateStateDependingOnManipFrame()
{
    if (this->useManipFrame)
        if (this->state != HIDDEN)
            this->setState(this->isManipulated() ? ACTIVE : INACTIVE);
}

void ControlPoint::checkIfGrabsMouse(int x, int y, const qglviewer::Camera * const cam)
{
    if (this->isManipulated())
    {
        setGrabsMouse(true);
        return;
    }
    if (this->state == HIDDEN || this->mesh.isHidden()) {
        setGrabsMouse(false);
        return;
    }

    qglviewer::Vec orig, dir;
    cam->convertClickToLine(QPoint(x,y), orig, dir);
    Vector3 rayOrigin = Vector3(orig), rayDir = Vector3(dir);

    auto [rotationAxis, intersectionPoint] = this->mouseOnRotationCircle(rayOrigin, rayDir);

    if (this->mouseOnCentralSphere(rayOrigin, rayDir)) {
        setGrabsMouse(true);
        this->currentAxis = NONE;
        this->isApplyingFreeMove = true;
        this->isApplyingRotation = false;
        this->isApplyingTranslation = false;
    }
    else if (this->mouseOnTranslationArrow(rayOrigin, rayDir)) {
        setGrabsMouse(true);
        this->isApplyingFreeMove = false;
        this->isApplyingRotation = false;
        this->isApplyingTranslation = true;
    } else if (rotationAxis != NONE) {
        setGrabsMouse(true);
        this->isApplyingFreeMove = false;
        this->isApplyingRotation = true;
        this->isApplyingTranslation = false;
        this->currentMousePosOnAction = intersectionPoint;
        this->currentAxis = rotationAxis;
    } else {
        setGrabsMouse(false);
        this->isApplyingFreeMove = false;
        this->isApplyingRotation = false;
        this->isApplyingTranslation = false;
    }


    // Constraints :
    qglviewer::WorldConstraint* constraint = new qglviewer::WorldConstraint();
    if (this->isApplyingFreeMove) {
        constraint->setTranslationConstraintType(qglviewer::AxisPlaneConstraint::FREE);
        constraint->setRotationConstraintType(qglviewer::AxisPlaneConstraint::FREE);
    } else if (this->isApplyingRotation) {
        constraint->setRotationConstraintType(qglviewer::AxisPlaneConstraint::AXIS);
        if (this->currentAxis == X)
            constraint->setRotationConstraintDirection(Vector3(1.0, 0.0, 0.0));
        else if (this->currentAxis == Y)
            constraint->setRotationConstraintDirection(Vector3(0.0, 1.0, 0.0));
        else if (this->currentAxis == Z)
            constraint->setRotationConstraintDirection(Vector3(0.0, 0.0, 1.0));
    } else if (this->isApplyingTranslation) {
        constraint->setTranslationConstraintType(qglviewer::AxisPlaneConstraint::AXIS);
        if (this->currentAxis == X)
            constraint->setTranslationConstraintDirection(Vector3(1.0, 0.0, 0.0));
        else if (this->currentAxis == Y)
            constraint->setTranslationConstraintDirection(Vector3(0.0, 1.0, 0.0));
        else if (this->currentAxis == Z)
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
    this->shape.position = this->getPosition();
    this->shape.radius = this->radius;
    this->shape.buildVerticesFlat();
    this->mesh.fromArray(this->shape.mesh.vertexArrayFloat);
    this->mesh.update();

    this->arrowSize = 3 * this->radius;
    this->circleRadius = 2 * this->radius;
}

void ControlPoint::display()
{
    if (this->state == HIDDEN || this->stillOnInitialState) return;
    else {
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

        if (this->mesh.shader != nullptr)
            this->mesh.shader->setVector("color", ControlPoint::GrabberStateColor[this->state]);

        if (this->translationMeshes.shader != nullptr) {
            // Display X (red)
            if (this->allowedTranslations[X]) {
                this->translationMeshes.shader->setVector("color", std::vector<float>({1.0, 0.0, 0.0, 1.0}));
                this->translationMeshes.fromArray({this->getPosition() - Vector3(1.0, 0.0, 0.0) * arrowSize, this->getPosition() + Vector3(1.0, 0.0, 0.0) * arrowSize});
                this->translationMeshes.display(GL_LINES, (isApplyingTranslation && currentAxis == X ? controlAxisSizeSelected : controlAxisSizeUnselected));
            }
            // Display Y (green)
            if (this->allowedTranslations[X]) {
                this->translationMeshes.shader->setVector("color", std::vector<float>({0.0, 1.0, 0.0, 1.0}));
                this->translationMeshes.fromArray({this->getPosition() - Vector3(0.0, 1.0, 0.0) * arrowSize, this->getPosition() + Vector3(0.0, 1.0, 0.0) * arrowSize});
                this->translationMeshes.display(GL_LINES, (isApplyingTranslation && currentAxis == Y ? controlAxisSizeSelected : controlAxisSizeUnselected));
            }
            // Display Z (blue)
            if (this->allowedTranslations[X]) {
                this->translationMeshes.shader->setVector("color", std::vector<float>({0.0, 0.0, 1.0, 1.0}));
                this->translationMeshes.fromArray({this->getPosition() - Vector3(0.0, 0.0, 1.0) * arrowSize, this->getPosition() + Vector3(0.0, 0.0, 1.0) * arrowSize});
                this->translationMeshes.display(GL_LINES, (isApplyingTranslation && currentAxis == Z ? controlAxisSizeSelected : controlAxisSizeUnselected));
            }
        } else if (this->mesh.shader != nullptr ){
            this->translationMeshes.shader = std::make_shared<Shader>(*this->mesh.shader);
        }
        if (this->rotationMeshes.shader != nullptr) {
            // Display X (red)
            if (this->allowedRotations[X]) {
                this->rotationMeshes.shader->setVector("color", std::vector<float>({1.0, 0.0, 0.0, 1.0}));
                this->rotationMeshes.fromArray(computeCircle(X));
                this->rotationMeshes.display(GL_LINES, (isApplyingRotation && currentAxis == X ? controlAxisSizeSelected : controlAxisSizeUnselected));
                if (isApplyingRotation && currentAxis == X) {
//                    this->rotationHelperSphere.translate(Vector3(1, 0, 0));
                    this->rotationHelperSphere.display();
                }
            }
            // Display Y (green)
            if (this->allowedRotations[Y]) {
                this->rotationMeshes.shader->setVector("color", std::vector<float>({0.0, 1.0, 0.0, 1.0}));
                this->rotationMeshes.fromArray(computeCircle(Y));
                this->rotationMeshes.display(GL_LINES, (isApplyingRotation && currentAxis == Y ? controlAxisSizeSelected : controlAxisSizeUnselected));
                if (isApplyingRotation && currentAxis == Y) {
//                    this->rotationHelperSphere.translate(Vector3(1, 0, 0));
                    this->rotationHelperSphere.display();
                }
            }
            // Display Z (blue)
            if (this->allowedRotations[Z]) {
                this->rotationMeshes.shader->setVector("color", std::vector<float>({0.0, 0.0, 1.0, 1.0}));
                this->rotationMeshes.fromArray(computeCircle(Z));
                this->rotationMeshes.display(GL_LINES, (isApplyingRotation && currentAxis == Z ? controlAxisSizeSelected : controlAxisSizeUnselected));
                if (isApplyingRotation && currentAxis == Z) {
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

Vector3 ControlPoint::getRotation() const
{
    if (!this->pressedPosBeforeAction.isValid())
        return Vector3(0.f, 0.f, 0.f);
    Vector3 rotation = this->pressedPosBeforeAction.getAllAnglesWith(this->currentMousePosOnAction);
    return rotation;
}

void ControlPoint::mousePressEvent(QMouseEvent * const event, qglviewer::Camera * const cam)
{
    this->pressedPosBeforeAction.setValid(false);
    if (this->grabsMouse()) {
        if (this->isApplyingRotation) {
            this->pressedPosBeforeAction = currentMousePosOnAction;
            Sphere s(1.f, this->getPosition() + circleRadius * this->pressedPosBeforeAction.normalized(), 6, 6);
            s.buildVerticesFlat();
            this->rotationHelperSphere.fromArray(s.mesh.vertexArray);
            this->startAction(QGLViewer::ROTATE);
        } else if (this->isApplyingTranslation) {
            this->startAction(QGLViewer::TRANSLATE);
        } else if (this->isApplyingFreeMove) {
            this->startAction(QGLViewer::TRANSLATE); // force translation
        }
    }
    qglviewer::ManipulatedFrame::mousePressEvent(event, cam);
}

void ControlPoint::mouseReleaseEvent(QMouseEvent * const event, qglviewer::Camera * const cam)
{
    qglviewer::ManipulatedFrame::mouseReleaseEvent(event, cam);
}

void ControlPoint::mouseMoveEvent(QMouseEvent * const event, qglviewer::Camera * const cam)
{
    QPoint mousePos = event->pos();
    qglviewer::Vec orig, dir;
    cam->convertClickToLine(mousePos, orig, dir);
    Vector3 rayOrigin = Vector3(orig), rayDir = Vector3(dir);
    this->currentMousePosOnAction = this->getIntersectionWithPlane(rayOrigin, rayDir, this->currentAxis) - this->getPosition();

    if (this->isApplyingRotation) {
        Sphere s(1.f, this->getPosition() + circleRadius * this->currentMousePosOnAction.normalized(), 6, 6);
        s.buildVerticesFlat();
        this->rotationHelperSphere.fromArray(s.mesh.vertexArray);
    }
    event->accept();
    qglviewer::ManipulatedFrame::mouseMoveEvent(event, cam);
}

void ControlPoint::wheelEvent(QWheelEvent * const event, qglviewer::Camera * const camera)
{
    setSphereRadius(this->radius - event->angleDelta().y()/10.f);
    this->startAction(QGLViewer::MouseAction::NO_MOUSE_ACTION);
    qglviewer::ManipulatedFrame::wheelEvent(event, camera);
}

void ControlPoint::setSphereRadius(float newRadius)
{
    if (minSphereRadius >= 0)
        newRadius = std::max(minSphereRadius, newRadius);
    if (maxSphereRadius >= 0)
        newRadius = std::min(maxSphereRadius, newRadius);
    this->radius = newRadius;
    this->updateSphere();
}

void ControlPoint::allowAllAxisTranslation(bool allow)
{
    this->allowedTranslations[X] = allow;
    this->allowedTranslations[Y] = allow;
    this->allowedTranslations[Z] = allow;
}

void ControlPoint::allowAllAxisRotations(bool allow)
{
    this->allowedRotations[X] = allow;
    this->allowedRotations[Y] = allow;
    this->allowedRotations[Z] = allow;
}

std::vector<Vector3> ControlPoint::computeCircle(Axis axis)
{
    std::vector<Vector3> points;
    for (int i = 0; i <= 360; i += 5) {
        float angle = i * PI / 180.f;
        float nextAngle = (i + 5) * PI / 180.f;
        if (axis == X) {
            points.push_back(this->getPosition() + Vector3(0.0, 1.0, 0.0).rotate(angle, 0, 0) * this->circleRadius);
            points.push_back(this->getPosition() + Vector3(0.0, 1.0, 0.0).rotate(nextAngle, 0, 0) * this->circleRadius);
        } else if (axis == Y) {
            points.push_back(this->getPosition() + Vector3(1.0, 0.0, 0.0).rotate(0, angle, 0) * this->circleRadius);
            points.push_back(this->getPosition() + Vector3(1.0, 0.0, 0.0).rotate(0, nextAngle, 0) * this->circleRadius);
        } else if (axis == Z) {
            points.push_back(this->getPosition() + Vector3(1.0, 0.0, 0.0).rotate(0, 0, angle) * this->circleRadius);
            points.push_back(this->getPosition() + Vector3(1.0, 0.0, 0.0).rotate(0, 0, nextAngle) * this->circleRadius);
        }
    }
    return points;
}

Vector3 ControlPoint::getIntersectionWithPlane(const Vector3& rayOrigin, const Vector3& rayDir, Axis axis)
{
    return Collision::intersectionRayPlane(rayOrigin, rayDir, this->getPosition(), Vector3(
                                               (axis == X ? 1.f :  0.f),
                                               (axis == Y ? 1.f :  0.f),
                                               (axis == Z ? 1.f :  0.f)));
}

bool ControlPoint::mouseOnCentralSphere(const Vector3& rayOrigin, const Vector3& rayDir)
{
    return Collision::intersectionRaySphere(rayOrigin, rayDir, this->getPosition(), this->radius).isValid();
}

bool ControlPoint::mouseOnTranslationArrow(const Vector3& rayOrigin, const Vector3& rayDir)
{
    float tolerence = this->radius/5.f;
    // X-axis
    if (this->allowedTranslations[X] && Collision::shortestDistanceBetweenSegments(rayOrigin, rayOrigin + rayDir * 1000.f,
                                        this->getPosition() - Vector3(1.0, 0.0, 0.0) * arrowSize,
                                        this->getPosition() + Vector3(1.0, 0.0, 0.0) * arrowSize) < tolerence) {
        this->currentAxis = X;
        return true;
    }
    // Y-axis
    if (this->allowedTranslations[Y] && Collision::shortestDistanceBetweenSegments(rayOrigin, rayOrigin + rayDir * 1000.f,
                                        this->getPosition() - Vector3(0.0, 1.0, 0.0) * arrowSize,
                                        this->getPosition() + Vector3(0.0, 1.0, 0.0) * arrowSize) < tolerence) {
        this->currentAxis = Y;
        return true;
    }
    // Z-axis
    if (this->allowedTranslations[Z] && Collision::shortestDistanceBetweenSegments(rayOrigin, rayOrigin + rayDir * 1000.f,
                                        this->getPosition() - Vector3(0.0, 0.0, 1.0) * arrowSize,
                                        this->getPosition() + Vector3(0.0, 0.0, 1.0) * arrowSize) < tolerence) {
        this->currentAxis = Z;
        return true;
    }
    return false;
}

std::pair<ControlPoint::Axis, Vector3> ControlPoint::mouseOnRotationCircle(const Vector3& rayOrigin, const Vector3& rayDir)
{
//    float circleRadiusSq = this->circleRadius * this->circleRadius;
    float tolerence = this->radius / 5.f;
    float minCircleRadiusSq = (this->circleRadius - tolerence) * (this->circleRadius - tolerence);
    float maxCircleRadiusSq = (this->circleRadius + tolerence) * (this->circleRadius + tolerence);
    Vector3 intersection;
    Vector3 bestIntersection;
    float distanceToCamSq = std::numeric_limits<float>::max();
    float distanceToCenterSq;
    Axis bestAxis = NONE;

    // X-axis
    if (this->allowedRotations[X]) {
        intersection = this->getIntersectionWithPlane(rayOrigin, rayDir, X);
        if (intersection.isValid()) {
            distanceToCenterSq = (intersection - this->getPosition()).norm2();
            if (minCircleRadiusSq < distanceToCenterSq && distanceToCenterSq < maxCircleRadiusSq) {
                if ((intersection - rayOrigin).norm2() < distanceToCamSq) {
                    distanceToCamSq = (intersection - rayOrigin).norm2();
                    bestAxis = X;
                    bestIntersection = intersection;
                }
            }
        }
    }
    // Y-axis
    if (this->allowedRotations[Y]) {
        intersection = this->getIntersectionWithPlane(rayOrigin, rayDir, Y);
        if (intersection.isValid()) {
            distanceToCenterSq = (intersection - this->getPosition()).norm2();
            if (minCircleRadiusSq < distanceToCenterSq && distanceToCenterSq < maxCircleRadiusSq) {
                if ((intersection - rayOrigin).norm2() < distanceToCamSq) {
                    distanceToCamSq = (intersection - rayOrigin).norm2();
                    bestAxis = Y;
                    bestIntersection = intersection;
                }
            }
        }
    }
    // Z-axis
    if (this->allowedRotations[Z]) {
        intersection = this->getIntersectionWithPlane(rayOrigin, rayDir, Z);
        if (intersection.isValid()) {
            distanceToCenterSq = (intersection - this->getPosition()).norm2();
            if (minCircleRadiusSq < distanceToCenterSq && distanceToCenterSq < maxCircleRadiusSq) {
                if ((intersection - rayOrigin).norm2() < distanceToCamSq) {
    //                distanceToCamSq = (intersection - rayOrigin).norm2();
                    bestAxis = Z;
                    bestIntersection = intersection;
                }
            }
        }
    }
/*
    if (bestAxis == NONE) {
        return false;
    } else {
        this->pressedPosBeforeAction = bestIntersection - this->getPosition();
        this->currentAxis = bestAxis;
        return true;
    }*/
    return {bestAxis, bestIntersection - this->getPosition()};
}

void ControlPoint::move(const Vector3& newPos)
{
    if (!this->isManipulated()) {
        this->setPosition(newPos.x, newPos.y, newPos.z);
        this->updateSphere();
    }
    this->stillOnInitialState = false;
}





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
