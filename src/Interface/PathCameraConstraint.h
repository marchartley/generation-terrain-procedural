#ifndef PATHCAMERACONSTRAINT_H
#define PATHCAMERACONSTRAINT_H

#include "Utils/Globals.h"
#include <QGLViewer/constraint.h>
#include <QGLViewer/camera.h>
#include "Utils/BSpline.h"
#include <vector>

class PathCameraConstraint : public qglviewer::Constraint
{
public:
    PathCameraConstraint(qglviewer::Camera* camera);
    PathCameraConstraint(qglviewer::Camera* camera, BSpline path);
    PathCameraConstraint(qglviewer::Camera* camera, std::vector<BSpline> paths);

    virtual void constrainTranslation(qglviewer::Vec& t, qglviewer::Frame* const fr);
    virtual void constrainRotation(qglviewer::Quaternion& q, qglviewer::Frame* const fr);

private:
    qglviewer::AxisPlaneConstraint* constraint;
    qglviewer::AxisPlaneConstraint* camRotationConstraint;
    qglviewer::Camera* camera;
    std::vector<BSpline> paths;
};

#endif // PATHCAMERACONSTRAINT_H
