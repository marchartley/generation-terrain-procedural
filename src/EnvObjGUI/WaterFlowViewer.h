#ifndef WATERFLOWVIEWER_H
#define WATERFLOWVIEWER_H

#include "GUIElements/ImageViewer.h"
#include "GUIElements/PainterToolsUI.h"

class WaterFlowViewer : public ImageViewer
{
    Q_OBJECT
protected: // Singleton
    WaterFlowViewer(const std::string& name, QWidget* parent = nullptr);
    WaterFlowViewer(const std::string& name, ChartView* chartView, QWidget* parent = nullptr);

public:
    DECLARE_PLOTTER_GETTER(WaterFlowViewer)

    virtual WaterFlowViewer* updateToolsInterface();

    // PainterToolParams painterParams;
    KelvinletToolParams kelvinletParams;

Q_SIGNALS:
    // void imagePainted(const GridF& newImage);
};

#endif // WATERFLOWVIEWER_H
