#include "FocusAreaViewer.h"

#include "GUIElements/PainterToolsUI.h"
#include "GUIElements/ImageViewerOptionsUI.h"

FocusAreaViewer::FocusAreaViewer(const std::string& name, QWidget* parent) : FocusAreaViewer(name, new ChartView(new Chart()), parent)
{}
FocusAreaViewer::FocusAreaViewer(const std::string& name, ChartView* chartView, QWidget* parent)
    : ImageViewer(name, chartView, parent)
{
    painterParams.additiveMode = true;
    painterParams.RGBimage = false;
    painterParams.minClampColor = Vector3(0, 0, 0);

    dataModel->imageData.displayParameters.colorRamp = BSpline({Vector3::red, Vector3::white, Vector3::green});

    QObject::connect(this, &FocusAreaViewer::movedOnImage, this, [&](const Vector3& clickPos, const Vector3& _prevPos, QMouseEvent* event) {
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
/*
FocusAreaViewer *FocusAreaViewer::getInstance(std::string name)
{
    if (name == "") name = FocusAreaViewer::defaultName;
    if (FocusAreaViewer::instances.count(name) == 0) {
        FocusAreaViewer::init(name);
    }
    return dynamic_cast<FocusAreaViewer*>(FocusAreaViewer::instances[name]);
}

FocusAreaViewer *FocusAreaViewer::init(const std::string& name, ChartView *chartView, QWidget *parent)
{
    if (FocusAreaViewer::instances.count(name))
        delete FocusAreaViewer::instances[name];
    FocusAreaViewer::instances[name] = new FocusAreaViewer(name, chartView, parent);
    return FocusAreaViewer::getInstance(name);
}
*/
FocusAreaViewer *FocusAreaViewer::updateToolsInterface()
{
    this->toolsInterface->clear();
    this->toolsInterface->add(PainterToolsUI::createPainterToolsUI(this->chartView, this->dataModel, &this->painterParams));
    ImageViewer::updateToolsInterface();
    return this;
}

/*
FocusAreaViewer *FocusAreaViewer::updateViewOptionsInterface()
{
    if (this->viewOptionsInterface != nullptr)
        this->viewOptionsInterface->clear();
    else
        this->viewOptionsInterface = new InterfaceUI(new QVBoxLayout());
    this->viewOptionsInterface->add(ImageViewerOptionsUI::createGreyImageViewerOptions(this->chartView, this->dataModel));
    return this;
}
*/
