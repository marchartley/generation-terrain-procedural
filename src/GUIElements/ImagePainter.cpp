#include "ImagePainter.h"

#include "GUIElements/PainterToolsUI.h"
#include "GUIElements/ImageViewerOptionsUI.h"

ImagePainter::ImagePainter(const std::string& name, QWidget* parent)
    : ImageViewer(name, parent)
{
    painterParams.additiveMode = true;
    painterParams.RGBimage = true;
    // painterParams.minClampColor = Vector3(0, 0, 0);

    // dataModel->imageData.displayParameters.colorRamp = BSpline({Vector3::white, Vector3::green});

    setOnMouseMoved([=](const Vector3& clickPos, const Vector3& _prevPos, QMouseEvent* event) {
        bool leftPressed = event->buttons().testFlag(Qt::LeftButton);
        bool rightPressed = event->buttons().testFlag(Qt::RightButton);
        if (!leftPressed && !rightPressed) return;
        if (this->painterParams.RGBimage) {
            GridV3 img = this->dataModel->getImage();
            PainterToolsUI::paintImage(img, clickPos, painterParams, rightPressed);
            Q_EMIT this->imagePainted(img);
            this->addImage(img);
        } else {
            GridF img = this->dataModel->getImageGrey();
            PainterToolsUI::paintImage(img, clickPos, painterParams, rightPressed);
            Q_EMIT this->imagePainted(img);
            this->addImage(img);
        }
        this->show();
    });
}

ImagePainter& ImagePainter::updateToolsInterface()
{
    this->toolsInterface->clear();
    auto UI = PainterToolsUI::createPainterToolsUI(&this->painterParams);
    this->toolsInterface->add(std::move(UI));
    ImageViewer::updateToolsInterface();
    return *this;
}

ImagePainter& ImagePainter::updateViewOptionsInterface()
{
    this->viewOptionsInterface->clear();

    if (this->painterParams.RGBimage)
        this->viewOptionsInterface->add(std::move(ImageViewerOptionsUI::createRGBImageViewerOptions(this)));
    else
        this->viewOptionsInterface->add(std::move(ImageViewerOptionsUI::createGreyImageViewerOptions(this)));

    return *this;
}
