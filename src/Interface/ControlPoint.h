#ifndef CONTROLPOINT_H
#define CONTROLPOINT_H

#include "Graphics/Sphere.h"
#include "Graphics/Mesh.h"
#include "DataStructure/Vector3.h"

#include <QGLViewer/manipulatedFrame.h>

enum GrabberState {
    HIDDEN   = 0x0,
    INACTIVE = 0x1,
    ACTIVE   = 0x2
};

class ControlPoint : public Mesh
{
public:
    ControlPoint();
    ControlPoint(Vector3 pos, float radius = 1.f, GrabberState state = INACTIVE, bool useTheManipulatedFrame = true);
    ~ControlPoint();

    void setState(GrabberState newState);

    void updateSphere();
    void move(Vector3 newPos);
    void display();

    Vector3 position;
    float radius;
    GrabberState state;
    bool useManipFrame;

    Sphere shape;

    qglviewer::ManipulatedFrame manipFrame;

    static std::shared_ptr<Shader> base_shader;
    static std::map<GrabberState, std::vector<float>> GrabberStateColor;
};

#endif // CONTROLPOINT_H
