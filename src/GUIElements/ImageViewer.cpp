#include "ImageViewer.h"
#include "GUIElements/ImageViewerOptionsUI.h"

ImageViewer::ImageViewer(const std::string& name, QWidget *parent) : ImageViewer(name, new ChartView(new Chart()), parent)
{
}

ImageViewer::ImageViewer(const std::string& name, ChartView *chartView, QWidget *parent) : AbstractPlotter(name, chartView, name, parent)
{
}
/*
ImageViewer *ImageViewer::getInstance(std::string name)
{
    if (name == "") name = ImageViewer::defaultName;
    if (ImageViewer::instances.count(ImageViewer::getIDname<ImageViewer>(name)) == 0) {
        ImageViewer::init(name);
    }
    return dynamic_cast<ImageViewer*>(ImageViewer::instances[ImageViewer::getIDname<ImageViewer>(name)]);
}

void ImageViewer::init(const std::string& name, ChartView *chartView, QWidget *parent)
{
    if (ImageViewer::instances.count(ImageViewer::getIDname<ImageViewer>(name)))
        delete ImageViewer::instances[ImageViewer::getIDname<ImageViewer>(name)];
    ImageViewer::instances[ImageViewer::getIDname<ImageViewer>(name)] = new ImageViewer(name, chartView, parent);
}
*/

ImageViewer* ImageViewer::updateUI()
{
    AbstractPlotter::updateUI();
    return this;
}

ImageViewer* ImageViewer::setNormalizedModeImage(bool normalize)
{
    this->dataModel->imageData.setNormalized(normalize);
    return updateUI();
}

ImageViewer* ImageViewer::setAbsoluteModeImage(bool absolute)
{
    this->dataModel->imageData.setAbsolute(absolute);
    return updateUI();
}

ImageViewer *ImageViewer::setFilteredValuesImage(bool filtered)
{
    this->dataModel->imageData.setClamped(filtered);
    return updateUI();
}

ImageViewer *ImageViewer::updateViewOptionsInterface()
{
    if (this->viewOptionsInterface != nullptr)
        this->viewOptionsInterface->clear();
    else
        this->viewOptionsInterface = new InterfaceUI(new QVBoxLayout());
    if (this->dataModel->imageData.image.isColor())
        this->viewOptionsInterface->add(ImageViewerOptionsUI::createRGBImageViewerOptions(this->chartView, this->dataModel));
    else
        this->viewOptionsInterface->add(ImageViewerOptionsUI::createGreyImageViewerOptions(this->chartView, this->dataModel));
    return this;
}


ImageViewer* ImageViewer::displayInfoUnderMouse(const Vector3 &relativeMousePos)
{
    if (relativeMousePos.minComp() < 0.f || relativeMousePos.maxComp() > 1.f)
        return this;
    std::ostringstream oss;
    if (this->hasImage()) {
        Vector3 size = this->dataModel->getImage().getDimensions();
        Vector3 position = relativeMousePos * size;
        Vector3 value = this->dataModel->getImage()[position];
        oss << "Image (" << int(position.x()) << ", " << int(position.y()) << ") = ";
        if (this->dataModel->imageData.image.isColor())
            oss << "(" << value.x() << ", " << value.y() << ", " << value.z() << ") ";
        else
            oss << value.x();
    }
    if (this->hasImage() && this->hasVectorField())
        oss << " -- ";
    if (this->hasVectorField()) {
        Vector3 size = this->dataModel->vectorData.field.getDimensions();
        Vector3 position = relativeMousePos * size;
        Vector3 value = this->dataModel->vectorData.getField()[position];
        oss << "Vector (" << int(position.x()) << ", " << int(position.y()) << ") = (" << value.x() << ", " << value.y() << ", " << value.z() << ") ";
    }
    this->mouseInfoLabel->setText(QString::fromStdString(oss.str()));
    return this;
}
