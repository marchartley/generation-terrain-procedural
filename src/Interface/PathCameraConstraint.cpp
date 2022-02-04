#include "PathCameraConstraint.h"

#include <QGLViewer/manipulatedCameraFrame.h>


PathCameraConstraint::PathCameraConstraint(qglviewer::Camera* camera)
    : PathCameraConstraint(camera, std::vector<BSpline>())
{
}

PathCameraConstraint::PathCameraConstraint(qglviewer::Camera* camera, BSpline path)
    : PathCameraConstraint(camera, std::vector<BSpline>(1, path))
{
}

PathCameraConstraint::PathCameraConstraint(qglviewer::Camera* camera, std::vector<BSpline> paths)
    : camera(camera), paths(paths)
{
    this->constraint = new qglviewer::WorldConstraint();
    this->constraint->setTranslationConstraintType(qglviewer::AxisPlaneConstraint::AXIS);
    this->camRotationConstraint = new qglviewer::CameraConstraint(camera);
    this->camRotationConstraint->setRotationConstraintType(qglviewer::AxisPlaneConstraint::PLANE);
    this->camRotationConstraint->setRotationConstraintDirection(Vector3(0, 0, 1));
//    this->constraint->setRotationConstraintType(qglviewer::AxisPlaneConstraint::FORBIDDEN);
}

void PathCameraConstraint::constrainTranslation(qglviewer::Vec &t, qglviewer::Frame * const fr)
{
    Vector3 vecT = Vector3(t.x, t.y, t.z);
    if (vecT.norm2() == 0)
        return; // Just a rotation, don't bother with translation computations
    std::vector<BSpline> closestPathCandidates;
    std::vector<int> closestPointIndexCandidates;
    BSpline closestPath;
    int closestPointIndex = -1;
    float dist = std::numeric_limits<float>::max();
    Vector3 currentPos = Vector3(fr->position().x, fr->position().y, fr->position().z);
    for (const BSpline& path : this->paths)
    {
        for (size_t i = 0; i < path.points.size(); i++)
        {
            float distToPoint = (path.points[i] - currentPos).norm2();
            if (std::abs(distToPoint - dist) < 1.0) {
                closestPathCandidates.push_back(path);
                closestPointIndexCandidates.push_back(i);
            } else if (distToPoint < dist) {
                dist = distToPoint;
                closestPathCandidates = std::vector<BSpline>(1, path);
                closestPointIndexCandidates = std::vector<int>(1, i);
            }
        }
    }
    bool goingForward = true;
    if (closestPathCandidates.size() == 1) {
        // If only one closest point, no doubt to have
        closestPath = closestPathCandidates[0];
        closestPointIndex = closestPointIndexCandidates[0];
    } else if (closestPathCandidates.size() > 1) {
        // Otherwise, choose the one for which the user is looking at
        float bestDotValue = -1.f;
        Vector3 camPos = this->camera->position();
        Vector3 camDir = this->camera->viewDirection();
        for (size_t i = 0; i < closestPathCandidates.size(); i++) {
            if ((closestPathCandidates[i].points[closestPointIndexCandidates[i]] - camPos).normalized().dot(camDir) > bestDotValue) {
                closestPath = closestPathCandidates[i];
                closestPointIndex = closestPointIndexCandidates[i];
                bestDotValue = (closestPath.points[closestPointIndex] - camPos).normalized().dot(camDir);
                goingForward = true;
            }
            if (closestPointIndexCandidates[i] > 0 &&
                    (closestPathCandidates[i].points[closestPointIndexCandidates[i] - 1] - camPos).normalized().dot(camDir) > bestDotValue) {
                closestPath = closestPathCandidates[i];
                closestPointIndex = closestPointIndexCandidates[i];
                bestDotValue = (closestPath.points[closestPointIndex - 1] - camPos).normalized().dot(camDir);
                goingForward = false;
            }
            if (closestPointIndexCandidates[i] < int(closestPathCandidates[i].points.size() - 1) &&
                    (closestPathCandidates[i].points[closestPointIndexCandidates[i] + 1] - camPos).normalized().dot(camDir) > bestDotValue) {
                closestPath = closestPathCandidates[i];
                closestPointIndex = closestPointIndexCandidates[i];
                bestDotValue = (closestPath.points[closestPointIndex + 1] - camPos).normalized().dot(camDir);
                goingForward = true;
            }
        }
    } else {
        return; // There is simply no path, we check all this for nothing...
    }
    if (!goingForward)
        t *= -1.f;

    float closestTimeOnCurve = closestPath.estimateClosestTime(this->camera->position());
    Vector3 constraintDir = closestPath.getDerivative(closestTimeOnCurve);
    this->constraint->setTranslationConstraintDirection(constraintDir);
    this->camera->setPosition(closestPath.getPoint(closestTimeOnCurve));
    this->camera->setPivotPoint(this->camera->position());
    this->camera->setViewDirection(closestPath.getDerivative(closestTimeOnCurve) * (goingForward || true ? 0.1f : -0.1f) + Vector3(this->camera->viewDirection()));
    this->camera->frame()->setSceneUpVector(Vector3(0, 0, 1));
    this->constraint->constrainTranslation(t, fr);
}

void PathCameraConstraint::constrainRotation(qglviewer::Quaternion &q, qglviewer::Frame * const fr)
{
//    this->camRotationConstraint-
    this->camRotationConstraint->constrainRotation(q, fr);
//    std::cout << "Rotation applied" << std::endl;
}
