#ifndef PAINTERTOOLSUI_H
#define PAINTERTOOLSUI_H

#include "GUIElements/CommonInterface.h"
#include "GUIElements/AbstractPlotter.h"

#include "DataStructure/Kelvinlet.h"

struct PainterToolParams {
    float radius = 10.f;
    float minRadius = 0.f;
    float maxRadius = 100.f;

    float opacity = 1.f;
    float minOpacity = 0.f;
    float maxOpacity = 1.f;

    bool additiveMode = true;

    float falloff = .5f;

    Vector3 color = Vector3(1, 1, 1);

    Vector3 minColor = Vector3(0, 0, 0);
    Vector3 maxColor = Vector3(1, 1, 1);

    Vector3 minClampColor = Vector3::min();
    Vector3 maxClampColor = Vector3::max();

    bool RGBimage = true;
};

struct KelvinletToolParams {
    std::vector<Kelvinlet*> kelvinlets;
    KelvinletPoint* currentKelvinlet = nullptr;

    float radialScale = 15.f; // epsilon
    float minRadialScale = 1.f;
    float maxRadialScale = 100.f;

    float mu = .5f; // 1.f; // Elastic shear modulus
    float minMu = .001f;
    float maxMu = 1.5f;

    float poisson = 0.5f; //.0f; // Poisson ratio (v = 1/2 is special case of Stokeslets [incompressible])
    float minPoisson = 0.f;
    float maxPoisson = 1.5f;

    float scale = 5.f; // For Scale Kelvinlet
    float minScale = 0.f;
    float maxScale = 100.f;

    Vector3 kelvinletPosition = Vector3::invalid;

    PlotVectorData temporaryVectorData;

    std::vector<Kelvinlet*> getKelvinlets() const { return this->kelvinlets; }
    GridV3 getInitialVectorField() const { return this->temporaryVectorData.field; }
    GridV3 getVectorField(bool takeIntoAccountCurrentKelvinlet = false) const;
};

class PainterToolsUI
{
public:
    static InterfaceUI* createPainterToolsUI(ChartView* chartView, PlotModel* dataModel, PainterToolParams* params);

    static GridV3& paintImage(GridV3& src, const Vector3& pos, PainterToolParams params, bool removeAmount = false);
    static GridF& paintImage(GridF& src, const Vector3& pos, PainterToolParams params, bool removeAmount = false);


    static InterfaceUI* createKelvinletToolsUI(ChartView* chartView, PlotModel* dataModel, KelvinletToolParams* params);
    static GridV3& paintKelvinlet(GridV3& src, const Vector3& pos, KelvinletToolParams* params);

protected:
    static void updateCurrentChartViewWithCurrentKelvinlets(ChartView *chartView, PlotModel *dataModel, KelvinletToolParams *params, const Vector3 &mouseRelPos, bool updateCurrentKelvinlet = true);
    static Kelvinlet* updateCurrentKelvinlet(KelvinletToolParams* params, const Vector3& mousePos);
};


#endif // PAINTERTOOLSUI_H
