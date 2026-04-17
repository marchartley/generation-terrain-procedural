#ifndef WATERFLOWVIEWER_H
#define WATERFLOWVIEWER_H

#include "GUIElements/ImageViewer.h"
#include "GUIElements/PainterToolsUI.h"

class WaterFlowViewer : public ImageViewer
{
    Q_OBJECT
public: // protected: // Singleton
    WaterFlowViewer(const std::string& name, QWidget* parent = nullptr);

public:
    DECLARE_PLOTTER_GETTER(WaterFlowViewer)

    virtual WaterFlowViewer& updateToolsInterface();

    // PainterToolParams painterParams;
    KelvinletToolParams kelvinletParams;

    DECLARE_EVENT(OnVectorFieldModified, (), ())
// Q_SIGNALS:
    // void vectorFieldModified();
    // void imagePainted(const GridF& newImage);
};

#endif // WATERFLOWVIEWER_H
