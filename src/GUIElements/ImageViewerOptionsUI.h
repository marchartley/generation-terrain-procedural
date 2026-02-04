#ifndef IMAGEVIEWEROPTIONSUI_H
#define IMAGEVIEWEROPTIONSUI_H

#include "GUIElements/CommonInterface.h"
#include "GUIElements/AbstractPlotter.h"

class ImageViewerOptionsUI
{
public:
    // ImageViewerOptionsUI();
    static InterfaceUI* createRGBImageViewerOptions(ChartView* chartView, PlotModel* dataModel);
    static InterfaceUI* createGreyImageViewerOptions(ChartView* chartView, PlotModel* dataModel);
};

#endif // IMAGEVIEWEROPTIONSUI_H
