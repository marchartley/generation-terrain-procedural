#include "EnvMaterialViewer.h"

#include "GUIElements/PainterToolsUI.h"
#include "GUIElements/ImageViewerOptionsUI.h"

EnvMaterialViewer::EnvMaterialViewer(const std::string& name, QWidget* parent) : EnvMaterialViewer(name, new ChartView(new Chart()), parent)
{}
EnvMaterialViewer::EnvMaterialViewer(const std::string& name, ChartView* chartView, QWidget* parent)
    : ImageViewer(name, chartView, parent)
{
    painterParams.additiveMode = true;
    painterParams.RGBimage = false;
    painterParams.minClampColor = Vector3(0, 0, 0);

    dataModel->imageData.displayParameters.colorRamp = BSpline({Vector3::white, Vector3::green});

    QObject::connect(this, &EnvMaterialViewer::movedOnImage, this, [&](const Vector3& clickPos, const Vector3& _prevPos, QMouseEvent* event) {
        bool leftPressed = event->buttons().testFlag(Qt::LeftButton);
        bool rightPressed = event->buttons().testFlag(Qt::RightButton);
        if (!leftPressed && !rightPressed) return;
        GridF img = this->dataModel->getImageGrey();
        PainterToolsUI::paintImage(img, clickPos, painterParams, rightPressed);
        Q_EMIT this->imagePainted(img);
        this->addImage(img);
        this->show();
    });
}

EnvMaterialViewer *EnvMaterialViewer::updateToolsInterface()
{
    this->toolsInterface->clear();
    this->toolsInterface->add(PainterToolsUI::createPainterToolsUI(this->chartView, this->dataModel, &this->painterParams));
    ImageViewer::updateToolsInterface();
    return this;
}
