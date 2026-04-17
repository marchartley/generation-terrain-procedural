#include "FocusAreaViewer.h"

#include "GUIElements/PainterToolsUI.h"
#include "GUIElements/ImageViewerOptionsUI.h"

FocusAreaViewer::FocusAreaViewer(const std::string& name, QWidget* parent)
    : ImageViewer(name, parent)
{
    painterParams.additiveMode = true;
    painterParams.RGBimage = false;
    painterParams.minClampColor = Vector3(0, 0, 0);

    dataModel->imageData.displayParameters.colorRamp = BSpline({Vector3::red, Vector3::white, Vector3::green});

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

FocusAreaViewer& FocusAreaViewer::updateToolsInterface()
{
    this->toolsInterface->clear();
    auto UI = PainterToolsUI::createPainterToolsUI(&this->painterParams);
    this->toolsInterface->add(std::move(UI));
    ImageViewer::updateToolsInterface();
    return *this;
}
