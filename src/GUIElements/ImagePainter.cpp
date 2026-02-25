#include "ImagePainter.h"

#include "GUIElements/PainterToolsUI.h"
#include "GUIElements/ImageViewerOptionsUI.h"

ImagePainter::ImagePainter(const std::string& name, QWidget* parent) : ImagePainter(name, new ChartView(new Chart()), parent)
{}
ImagePainter::ImagePainter(const std::string& name, ChartView* chartView, QWidget* parent)
    : ImageViewer(name, chartView, parent)
{
    painterParams.additiveMode = true;
    painterParams.RGBimage = true;
    // painterParams.minClampColor = Vector3(0, 0, 0);

    // dataModel->imageData.displayParameters.colorRamp = BSpline({Vector3::white, Vector3::green});

    QObject::connect(this, &ImagePainter::movedOnImage, this, [&](const Vector3& clickPos, const Vector3& _prevPos, QMouseEvent* event) {
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
/*
ImagePainter *ImagePainter::getInstance(const std::string& name)
{
    if (name == "") name = ImagePainter::defaultName;
    if (ImagePainter::instances.count(name) == 0) {
        ImagePainter::init(name);
    }
    return dynamic_cast<ImagePainter*>(ImagePainter::instances[name]);
}

ImagePainter *ImagePainter::init(const std::string& name, ChartView *chartView, QWidget *parent)
{
    if (ImagePainter::instances.count(name))
        delete ImagePainter::instances[name];
    ImagePainter::instances[name] = new ImagePainter(name, chartView, parent);
    return ImagePainter::getInstance(name);
}
*/
ImagePainter *ImagePainter::updateToolsInterface()
{
    this->toolsInterface->clear();
    this->toolsInterface->add(PainterToolsUI::createPainterToolsUI(this->chartView, this->dataModel, &this->painterParams));
    ImageViewer::updateToolsInterface();
    return this;
}

ImagePainter *ImagePainter::updateViewOptionsInterface()
{
    if (this->viewOptionsInterface != nullptr)
        this->viewOptionsInterface->clear();
    else
        this->viewOptionsInterface = new InterfaceUI(new QVBoxLayout());

    if (this->painterParams.RGBimage)
        this->viewOptionsInterface->add(ImageViewerOptionsUI::createRGBImageViewerOptions(this->chartView, this->dataModel));
    else
        this->viewOptionsInterface->add(ImageViewerOptionsUI::createGreyImageViewerOptions(this->chartView, this->dataModel));

    return this;
}
