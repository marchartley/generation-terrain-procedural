#include "EnvMaterialViewer.h"

#include "GUIElements/PainterToolsUI.h"
#include "GUIElements/ImageViewerOptionsUI.h"

EnvMaterialViewer::EnvMaterialViewer(const std::string& name, QWidget* parent)
    : ImageViewer(name, parent)
{
    painterParams.additiveMode = true;
    painterParams.RGBimage = false;
    painterParams.minClampColor = Vector3(0, 0, 0);

    dataModel->imageData.displayParameters.colorRamp = BSpline({Vector3::white, Vector3::green});

    setOnMouseMoved([=](const Vector3& clickPos, const Vector3& _prevPos, QMouseEvent* event) {
        bool leftPressed = event->buttons().testFlag(Qt::LeftButton);
        bool rightPressed = event->buttons().testFlag(Qt::RightButton);
        if (!leftPressed && !rightPressed) return;
        GridF img = this->dataModel->getImageGrey();
        PainterToolsUI::paintImage(img, clickPos, painterParams, rightPressed);
        // Q_EMIT this->imagePainted(img);
        emitOnImagePainted(img);
        this->addImage(img);
        this->show();
    });
}

EnvMaterialViewer& EnvMaterialViewer::updateToolsInterface()
{
    this->toolsInterface->clear();
    this->toolsInterface->add(PainterToolsUI::createPainterToolsUI(&this->painterParams));
    ImageViewer::updateToolsInterface();
    return *this;
}
