#ifndef VISITINGCAMERA_H
#define VISITINGCAMERA_H

#include "Utils/BSpline.h"
#include <QGLViewer/camera.h>

class VisitingCamera : public qglviewer::Camera
{
public:
    VisitingCamera();

    void moveForward(float distance);
    void moveBackward(float distance);

    bool isVisiting = false;

    std::vector<BSpline> paths;
};

#endif // VISITINGCAMERA_H
