#ifndef PAINTERTOOLSUI_H
#define PAINTERTOOLSUI_H

#include "GUIElements/CommonInterface.h"
#include "GUIElements/AbstractPlotter.h"

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

class PainterToolsUI
{
public:
    static InterfaceUI* createPainterToolsUI(ChartView* chartView, PlotModel* dataModel, PainterToolParams* params);
    static GridV3& paintImage(GridV3& src, const Vector3& pos, PainterToolParams params, bool removeAmount = false);
    static GridF& paintImage(GridF& src, const Vector3& pos, PainterToolParams params, bool removeAmount = false);
};


#endif // PAINTERTOOLSUI_H
