#include "ControlPoint.h"

#include "Utils/Utils.h"

std::shared_ptr<Shader> ControlPoint::base_shader = nullptr;
std::map<GrabberState, std::vector<float>> ControlPoint::default_GrabberStateColor = {
    {GrabberState::HIDDEN, {.0f, .0f, .0f, 0.f}},
    {GrabberState::INACTIVE, {.3f, .0f, .0f, .5f}},
    {GrabberState::ACTIVE, {.8f, .0f, .0f, .8f}},
    {GrabberState::POSITIVE, {.2f, 1.f, .1f, .8f}},
    {GrabberState::NEGATIVE, {1.f, .2f, .1f, 8.f}},
    {GrabberState::NEUTRAL, {.8f, .8f, .8f, .8f}},
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

ControlPoint::ControlPoint(Vector3 pos, float radius, GrabberState state, bool useTheManipulatedFrame)
    : /*CustomInteractiveObject(), pos(pos),*/ radius(radius), state(state), useManipFrame(useTheManipulatedFrame), shape(radius, getPosition(), 10, 10)
{
    this->mesh = Mesh((ControlPoint::base_shader ? std::make_shared<Shader>(*ControlPoint::base_shader) : nullptr), true);
    this->move(pos);
    this->prevPosition = pos;
    this->GrabberStateColor = ControlPoint::default_GrabberStateColor;
    if (!useManipFrame) {
        this->removeFromMouseGrabberPool();
    } else {
        this->addInMouseGrabberPool();
    }
    this->allowAllAxisRotations(false);
    this->allowAllAxisTranslation(false);
    this->rotationMeshes.show();
    this->translationMeshes.show();

    QObject::connect(this, &qglviewer::ManipulatedFrame::modified, this, [=](){
        Q_EMIT this->modified();
        if ((this->prevPosition - this->getPosition()).norm2() > 1.0) {
            this->prevPosition = this->getPosition();
//        if (this->positionsHistory.empty() || this->positionsHistory.back() != this->prevPosition) {
            this->positionsHistory.push_back(prevPosition);
            if (this->positionsHistory.size() > 10) {
                this->positionsHistory.erase(positionsHistory.begin(), std::max(positionsHistory.end() - 10, positionsHistory.begin()));
//                std::cout << positionsHistory.size() << " pos stored" << std::endl;
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
    } else if (this->mouseOnRotationCircle(rayOrigin, rayDir)) {
        setGrabsMouse(true);
        this->isApplyingFreeMove = false;
        this->isApplyingRotation = true;
        this->isApplyingTranslation = false;
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
    if(this->onUpdateCallback) { // && this->useManipFrame && this->/*manipFrame.*/isManipulated())
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
        Q_EMIT this->afterModified();
    }
    this->currentlyManipulated = this->isManipulated();

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
    if (this->state == HIDDEN) return;
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
            }
            // Display Y (green)
            if (this->allowedRotations[Y]) {
                this->rotationMeshes.shader->setVector("color", std::vector<float>({0.0, 1.0, 0.0, 1.0}));
                this->rotationMeshes.fromArray(computeCircle(Y));
                this->rotationMeshes.display(GL_LINES, (isApplyingRotation && currentAxis == Y ? controlAxisSizeSelected : controlAxisSizeUnselected));
            }
            // Display Z (blue)
            if (this->allowedRotations[Z]) {
                this->rotationMeshes.shader->setVector("color", std::vector<float>({0.0, 0.0, 1.0, 1.0}));
                this->rotationMeshes.fromArray(computeCircle(Z));
                this->rotationMeshes.display(GL_LINES, (isApplyingRotation && currentAxis == Z ? controlAxisSizeSelected : controlAxisSizeUnselected));
            }
        } else if (this->mesh.shader != nullptr ){
            this->rotationMeshes.shader = std::make_shared<Shader>(*this->mesh.shader);
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
//    CustomInteractiveObject::hide();
}

void ControlPoint::show()
{
    if (useManipFrame)
        this->addInMouseGrabberPool();
    this->mesh.show();
    this->translationMeshes.show();
    this->rotationMeshes.show();
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

void ControlPoint::mousePressEvent(QMouseEvent * const event, qglviewer::Camera * const cam)
{
    if (this->grabsMouse()) {
        if (this->isApplyingRotation) {
            this->startAction(QGLViewer::ROTATE);
        }
        else if (this->isApplyingTranslation) {
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
    qglviewer::ManipulatedFrame::mouseMoveEvent(event, cam);
}

void ControlPoint::wheelEvent(QWheelEvent * const event, qglviewer::Camera * const camera)
{
    setSphereRadius(this->radius - event->angleDelta().y()/10.f);
    this->startAction(QGLViewer::MouseAction::NO_MOUSE_ACTION);
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
    float pi = 3.141592;
    std::vector<Vector3> points;
    for (int i = 0; i <= 360; i += 5) {
        float angle = i * pi / 180.f;
        float nextAngle = (i + 5) * pi / 180.f;
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

bool ControlPoint::mouseOnCentralSphere(Vector3 rayOrigin, Vector3 rayDir)
{
    return intersectionRaySphere(rayOrigin, rayDir, this->getPosition(), this->radius).isValid();
}

bool ControlPoint::mouseOnTranslationArrow(Vector3 rayOrigin, Vector3 rayDir)
{
    float tolerence = this->radius/5.f;
    // X-axis
    if (this->allowedTranslations[X] && shortestDistanceBetweenSegments(rayOrigin, rayOrigin + rayDir * 1000.f,
                                        this->getPosition() - Vector3(1.0, 0.0, 0.0) * arrowSize,
                                        this->getPosition() + Vector3(1.0, 0.0, 0.0) * arrowSize) < tolerence) {
        this->currentAxis = X;
        return true;
    }
    // Y-axis
    if (this->allowedTranslations[Y] && shortestDistanceBetweenSegments(rayOrigin, rayOrigin + rayDir * 1000.f,
                                        this->getPosition() - Vector3(0.0, 1.0, 0.0) * arrowSize,
                                        this->getPosition() + Vector3(0.0, 1.0, 0.0) * arrowSize) < tolerence) {
        this->currentAxis = Y;
        return true;
    }
    // Z-axis
    if (this->allowedTranslations[Z] && shortestDistanceBetweenSegments(rayOrigin, rayOrigin + rayDir * 1000.f,
                                        this->getPosition() - Vector3(0.0, 0.0, 1.0) * arrowSize,
                                        this->getPosition() + Vector3(0.0, 0.0, 1.0) * arrowSize) < tolerence) {
        this->currentAxis = Z;
        return true;
    }
    return false;
}

bool ControlPoint::mouseOnRotationCircle(Vector3 rayOrigin, Vector3 rayDir)
{
//    float circleRadiusSq = this->circleRadius * this->circleRadius;
    float tolerence = this->radius / 5.f;
    float minCircleRadiusSq = (this->circleRadius - tolerence) * (this->circleRadius - tolerence);
    float maxCircleRadiusSq = (this->circleRadius + tolerence) * (this->circleRadius + tolerence);
    Vector3 intersection;
    float distanceToCamSq = std::numeric_limits<float>::max();
    float distanceToCenterSq;
    Axis bestAxis = NONE;

    // X-axis
    if (this->allowedRotations[X]) {
        intersection = intersectionRayPlane(rayOrigin, rayDir, this->getPosition(), Vector3(1.0, 0.0, 0.0));
        if (intersection.isValid()) {
            distanceToCenterSq = (intersection - this->getPosition()).norm2();
            if (minCircleRadiusSq < distanceToCenterSq && distanceToCenterSq < maxCircleRadiusSq) {
                if ((intersection - rayOrigin).norm2() < distanceToCamSq) {
                    distanceToCamSq = (intersection - rayOrigin).norm2();
                    bestAxis = X;
                }
            }
        }
    }
    // Y-axis
    if (this->allowedRotations[Y]) {
        intersection = intersectionRayPlane(rayOrigin, rayDir, this->getPosition(), Vector3(0.0, 1.0, 0.0));
        if (intersection.isValid()) {
            distanceToCenterSq = (intersection - this->getPosition()).norm2();
            if (minCircleRadiusSq < distanceToCenterSq && distanceToCenterSq < maxCircleRadiusSq) {
                if ((intersection - rayOrigin).norm2() < distanceToCamSq) {
                    distanceToCamSq = (intersection - rayOrigin).norm2();
                    bestAxis = Y;
                }
            }
        }
    }
    // Z-axis
    if (this->allowedRotations[Z]) {
        intersection = intersectionRayPlane(rayOrigin, rayDir, this->getPosition(), Vector3(0.0, 0.0, 1.0));
        if (intersection.isValid()) {
            distanceToCenterSq = (intersection - this->getPosition()).norm2();
            if (minCircleRadiusSq < distanceToCenterSq && distanceToCenterSq < maxCircleRadiusSq) {
                if ((intersection - rayOrigin).norm2() < distanceToCamSq) {
    //                distanceToCamSq = (intersection - rayOrigin).norm2();
                    bestAxis = Z;
                }
            }
        }
    }

    if (bestAxis == NONE) {
        return false;
    } else {
        this->currentAxis = bestAxis;
        return true;
    }
}

void ControlPoint::move(Vector3 newPos)
{
    this->setPosition(newPos.x, newPos.y, newPos.z);
    this->updateSphere();
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
