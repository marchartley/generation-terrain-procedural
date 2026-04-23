#include "WaterFlowViewer.h"

#include "GUIElements/ImageViewerOptionsUI.h"

WaterFlowViewer::WaterFlowViewer(const std::string& name, QWidget* parent)
    : ImageViewer(name, parent)
{
    dataModel->imageData.displayParameters.colorRamp = BSpline({Vector3::red, Vector3::white, Vector3::green});

    kelvinletParams.setOnNewKelvinlet([=](Kelvinlet*) { emitVectorFieldModified(); });
}

WaterFlowViewer& WaterFlowViewer::updateToolsInterface()
{
    this->toolsInterface->clear();
    auto UI = PainterToolsUI::createKelvinletToolsUI(this, &this->kelvinletParams);
    this->toolsInterface->add(std::move(UI));
    return *this;
}
