#include "ImageViewer.h"
#include "GUIElements/ImageViewerOptionsUI.h"

ImageViewer::ImageViewer(const std::string& name, QWidget *parent) : AbstractPlotter(name, name, parent)
{
}

ImageViewer& ImageViewer::updateUI()
{
    AbstractPlotter::updateUI();
    return *this;
}

ImageViewer& ImageViewer::setNormalizedModeImage(bool normalize)
{
    this->dataModel->imageData.setNormalized(normalize);
    return updateUI();
}

ImageViewer& ImageViewer::setAbsoluteModeImage(bool absolute)
{
    this->dataModel->imageData.setAbsolute(absolute);
    return updateUI();
}

ImageViewer& ImageViewer::setFilteredValuesImage(bool filtered)
{
    this->dataModel->imageData.setClamped(filtered);
    return updateUI();
}

ImageViewer& ImageViewer::updateViewOptionsInterface()
{
    this->viewOptionsInterface->clear();

    if (this->dataModel->imageData.image.isColor()) {
        this->viewOptionsInterface->add(std::move(ImageViewerOptionsUI::createRGBImageViewerOptions(this)));
    } else {
        this->viewOptionsInterface->add(std::move(ImageViewerOptionsUI::createGreyImageViewerOptions(this)));
    }
    return *this;
}


void ImageViewer::displayInfoUnderMouse(const Vector3 &relativeMousePos)
{
    if (relativeMousePos.minComp() < 0.f || relativeMousePos.maxComp() > 1.f)
        return;
    std::ostringstream oss;
    if (this->hasImage()) {
        Vector3 size = this->dataModel->getImage().getDimensions();
        Vector3 position = relativeMousePos * size;
        Vector3 value = this->dataModel->getImage()[position];
        oss << "Image (" << int(position.x()) << ", " << int(position.y()) << ") = ";
        if (this->dataModel->imageData.image.isColor())
            oss << "(" << value.x() << ", " << value.y() << ", " << value.z() << ")  [norm: " << value.norm() << "]";
        else
            oss << value.x();
    }
    if (this->hasImage() && this->hasVectorField())
        oss << " -- ";
    if (this->hasVectorField()) {
        Vector3 size = this->dataModel->vectorData.field.getDimensions();
        Vector3 position = relativeMousePos * size;
        Vector3 value = this->dataModel->vectorData.getField()[position];
        oss << "Vector (" << int(position.x()) << ", " << int(position.y()) << ") = (" << value.x() << ", " << value.y() << ", " << value.z() << ")  [norm: " << value.norm() << "]";
    }
    this->mouseInfoLabel->setText(QString::fromStdString(oss.str()));
}
