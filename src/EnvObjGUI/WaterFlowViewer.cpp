#include "WaterFlowViewer.h"

#include "GUIElements/ImageViewerOptionsUI.h"

WaterFlowViewer::WaterFlowViewer(const std::string& name, QWidget* parent) : WaterFlowViewer(name, new ChartView(new Chart()), parent)
{}
WaterFlowViewer::WaterFlowViewer(const std::string& name, ChartView* chartView, QWidget* parent)
    : ImageViewer(name, chartView, parent)
{
    dataModel->imageData.displayParameters.colorRamp = BSpline({Vector3::red, Vector3::white, Vector3::green});
}

WaterFlowViewer *WaterFlowViewer::updateToolsInterface()
{
    this->toolsInterface->clear();
    this->toolsInterface->add(PainterToolsUI::createKelvinletToolsUI(this, this->chartView, this->dataModel, &this->kelvinletParams));
    return this;
}
