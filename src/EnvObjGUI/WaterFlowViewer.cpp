#include "WaterFlowViewer.h"

#include "GUIElements/ImageViewerOptionsUI.h"

WaterFlowViewer::WaterFlowViewer(const std::string& name, QWidget* parent) : WaterFlowViewer(name, new ChartView(new Chart()), parent)
{}
WaterFlowViewer::WaterFlowViewer(const std::string& name, ChartView* chartView, QWidget* parent)
    : ImageViewer(name, chartView, parent)
{
    dataModel->imageData.displayParameters.colorRamp = BSpline({Vector3::red, Vector3::white, Vector3::green});

    QObject::connect(this, &WaterFlowViewer::movedOnImage, this, [&](const Vector3& clickPos, const Vector3& _prevPos, QMouseEvent* event) {
        bool leftPressed = event->buttons().testFlag(Qt::LeftButton);
        bool rightPressed = event->buttons().testFlag(Qt::RightButton);
        if (!leftPressed && !rightPressed) return;
        // GridF img = this->dataModel->getImageGrey();
        // PainterToolsUI::paintImage(img, clickPos, painterParams, rightPressed);
        // Q_EMIT this->imagePainted(img);
        // this->addImage(img);

        if (kelvinletParams.currentKelvinlet && kelvinletParams.kelvinletPosition.isValid()) {
            GridV3 newFlow = this->dataModel->vectorData.field;
            PainterToolsUI::paintKelvinlet(newFlow, clickPos, &this->kelvinletParams);
            auto [overlay, alpha] = PlotVectorData::createFieldImageAndAlpha(newFlow, Vector3i(200, 200, 1), Vector3i(20, 20, 1));
            this->setOverlay(overlay, alpha, "kelvinlet preview");
        }
        this->show();
    });
    QObject::connect(this, &WaterFlowViewer::clickedOnImage, this, [&](const Vector3& pos, const Vector3& value, bool leftClick, bool rightClick) {
        if (pos.isValid() && rightClick) { // Cancelled
            kelvinletParams.kelvinletPosition = Vector3::invalid();
            return;
        }
        if (!pos.isValid() && leftClick) { // Mouse released
            GridV3& newFlow = this->dataModel->vectorData.field;
            PainterToolsUI::paintKelvinlet(newFlow, pos, &this->kelvinletParams);
            this->hideOverlay("kelvinlet preview");
            this->show();
        }
    });
}

/*
WaterFlowViewer *WaterFlowViewer::getInstance(std::string name)
{
    if (name == "") name = WaterFlowViewer::defaultName;
    if (WaterFlowViewer::instances.count(name) == 0) {
        WaterFlowViewer::init(name);
    }
    return dynamic_cast<WaterFlowViewer*>(WaterFlowViewer::instances[name]);
}

WaterFlowViewer *WaterFlowViewer::init(const std::string& name, ChartView *chartView, QWidget *parent)
{
    if (WaterFlowViewer::instances.count(name))
        delete WaterFlowViewer::instances[name];
    WaterFlowViewer::instances[name] = new WaterFlowViewer(name, chartView, parent);
    return WaterFlowViewer::getInstance(name);
}
*/

WaterFlowViewer *WaterFlowViewer::updateToolsInterface()
{
    // this->toolsInterface->clear();
    // this->toolsInterface->add(PainterToolsUI::createPainterToolsUI(this->chartView, this->dataModel, &this->painterParams));
    // ImageViewer::updateToolsInterface();
    this->toolsInterface->clear();
    this->toolsInterface->add(PainterToolsUI::createKelvinletToolsUI(this->chartView, this->dataModel, &this->kelvinletParams));
    return this;
}
