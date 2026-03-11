#ifndef IMAGEVIEWEROPTIONSUI_H
#define IMAGEVIEWEROPTIONSUI_H

#include "GUIElements/CommonInterface.h"
#include "GUIElements/AbstractPlotter.h"

class ImageViewerOptionsUI
{
public:
    // ImageViewerOptionsUI();
    static InterfaceUI* createRGBImageViewerOptions(AbstractPlotter* plotter);
    static InterfaceUI* createGreyImageViewerOptions(AbstractPlotter* plotter);
};

#endif // IMAGEVIEWEROPTIONSUI_H
