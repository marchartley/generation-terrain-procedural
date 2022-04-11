#ifndef CONTROLPOINT_H
#define CONTROLPOINT_H

#include "Graphics/Sphere.h"
#include "Graphics/Mesh.h"
#include "DataStructure/Vector3.h"

#include "Interface/CustomInteractiveObject.h"
#include <QGLViewer/manipulatedFrame.h>

enum GrabberState {
    HIDDEN         = 0x0,
    INACTIVE       = 0x1,
    ACTIVE         = 0x2,
    POSITIVE       = 0x3,
    NEGATIVE       = 0x4,
    NEUTRAL        = 0x5,
    CUSTOM_STATE_0 = 0x6,
    CUSTOM_STATE_1 = 0x7,
    CUSTOM_STATE_2 = 0x8,
    CUSTOM_STATE_3 = 0x9,
    CUSTOM_STATE_4 = 0xA,
    CUSTOM_STATE_5 = 0xB,
    CUSTOM_STATE_6 = 0xC,
    CUSTOM_STATE_7 = 0xD,
    CUSTOM_STATE_8 = 0xE,
    CUSTOM_STATE_9 = 0xF,
};

class ControlPoint : public CustomInteractiveObject
{
    Q_OBJECT
public:
    ControlPoint();
    ControlPoint(Vector3 pos, float radius = 1.f, GrabberState state = INACTIVE, bool useTheManipulatedFrame = true);
    ~ControlPoint();

    void setState(GrabberState newState);
    void updateStateDependingOnManipFrame();

Q_SIGNALS:
    void modified();
    void afterModified();

public Q_SLOTS:
    void hide();
    void show();

private:
    void onUpdate(std::function<void()> func);
    void afterUpdate(std::function<void()> func);
public:
    void updateSphere();
    void move(Vector3 newPos);
    void display();

    void setGrabberStateColor(std::map<GrabberState, std::vector<float>> stateColorMap);
    void setGrabberStateColor(GrabberState state, std::vector<float> color);

    Vector3 position;
    Vector3 prevPosition;
    float radius;
    GrabberState state;
    bool useManipFrame;
    bool currentlyManipulated;

    Mesh mesh;

    Sphere shape;

    std::function<void()> onUpdateCallback;
    std::function<void()> afterUpdateCallback;

    qglviewer::ManipulatedFrame manipFrame;
    std::map<GrabberState, std::vector<float>> GrabberStateColor;

    static std::shared_ptr<Shader> base_shader;
    static std::map<GrabberState, std::vector<float>> default_GrabberStateColor;
};

#endif // CONTROLPOINT_H
