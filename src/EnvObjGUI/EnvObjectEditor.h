#ifndef ENVOBJECTEDITOR_H
#define ENVOBJECTEDITOR_H

#include "GUIElements/ImageViewer.h"
#include "GUIElements/PainterToolsUI.h"

#include "EnvObject/EnvObject.h"
#include "EnvObject/EnvPoint.h"
#include "EnvObject/EnvCurve.h"
#include "EnvObject/EnvArea.h"

struct BodyKelvinletParameters {
    float force = 1.f;
    float minForce = -100.f;
    float maxForce = 100.f;

    float radius = 5.f;
    float minRadius = 0.f;
    float maxRadius = 100.f;


    PinchKelvinletCurve* pinchK;
    TwistKelvinletCurve* twistK;
    GrabKelvinletCurve* grabK;
    ScaleKelvinletCurve* scaleK;

    std::vector<RelativeKelvinlet*> relativeKelvinlets;
};

class EnvObjectEditor : public ImageViewer
{
    Q_OBJECT
protected: // Singleton
    EnvObjectEditor(const std::string& name, QWidget* parent = nullptr);
    EnvObjectEditor(const std::string& name, ChartView* chartView, QWidget* parent = nullptr);

public:
    DECLARE_PLOTTER_GETTER(EnvObjectEditor)

    virtual EnvObjectEditor* updateToolsInterface();

    EnvObjectEditor* addEnvObject(EnvObject* envObj);

    std::pair<GridV3, GridF> displayEnvObject() const;

    GridF& simulateDeposition(GridF& currentState, int iterations = 10);

    GridV3 getVectorFieldWithRotation(bool takeIntoAccountCurrentKelvinlet) const;

    void updateCurrentChartViewWithCurrentKelvinlets(const Vector3& mouseRelPos, bool updateCurrentKelvinlet);

    EnvObject* validateEnvObject(bool takeIntoAccountCurrentKelvinlet = true) const;

    void animateEnvObject(bool animate);

    // PainterToolParams painterParams;
    KelvinletToolParams kelvinletParams;

    BodyKelvinletParameters bodyParameters;

    EnvObjectInstance* currentObject;
    GridF depositionGrid;
    bool depositionSimulationDisplay = false;

    enum KELVINLET_ANCHOR_POINT {
        MAIN, START, END, UNDEFINED
    };
    KELVINLET_ANCHOR_POINT currentAnchorPoint = UNDEFINED;
    std::map<Kelvinlet*, KELVINLET_ANCHOR_POINT> kelvinletAnchors;

    bool animating = false;
    std::vector<Vector3> verticesTargets;
    int animationFrame = 0;

    float objectScale = 1.f;

Q_SIGNALS:
           // void imagePainted(const GridF& newImage);
};

#endif // ENVOBJECTEDITOR_H
