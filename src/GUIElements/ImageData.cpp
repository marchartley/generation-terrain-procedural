#include "ImageData.h"


PlotImageData::PlotImageData() : PlotImageData(GridV3())
{}

PlotImageData::PlotImageData(const GridV3 &img) : image(img)
{}

PlotImageData::PlotImageData(const GridF &img) : image(img)
{}

PlotImageData& PlotImageData::setImage(const GridV3 &img)
{
    this->image.setImage(img);
    // this->callOnImageModifiedCallbacks();
    emitImageModified();
    return *this;
}

PlotImageData& PlotImageData::setImage(const GridF &img)
{
    this->image.setImage(img);
    // this->callOnImageModifiedCallbacks();
    emitImageModified();
    return *this;
}

PlotImageData& PlotImageData::setNormalized(bool normalize)
{
    this->displayParameters.normalized = normalize;
    return *this;
}

PlotImageData& PlotImageData::setColorRanges(const Vector3 &minRange, const Vector3 &maxRange)
{
    this->displayParameters.colorRangeMin = minRange;
    this->displayParameters.colorRangeMax = maxRange;
    return *this;
}

PlotImageData& PlotImageData::setAbsolute(bool absolute)
{
    this->displayParameters.absolute = absolute;
    return *this;
}

PlotImageData& PlotImageData::setClamped(bool clamp)
{
    this->displayParameters.clamped = clamp;
    return *this;
}

GridV3 PlotImageData::prepareImageForDisplay(const Image& img) const
{    if (img.empty()) return GridV3();
    auto displayedImage = img.getColorImage(); //.resize(imgSize);
    if (displayedImage.empty()) return GridV3();

    if (this->displayParameters.clamped) {
        displayedImage.iterateParallel([&](size_t i) {
            for (int c = 0; c < 3; c++) {
                displayedImage[i][c] = std::clamp(displayedImage[i][c], this->displayParameters.colorRangeMin[c], this->displayParameters.colorRangeMax[c]);
            }
        });
    }

    if (this->displayParameters.absolute) {
        displayedImage = displayedImage.abs();
    }

    if (this->displayParameters.normalized) {
        for (int c = 0; c < 3; c++) {
            float min = std::numeric_limits<float>::max();
            float max = std::numeric_limits<float>::lowest();
            displayedImage.iterate([&](size_t i) {
                min = std::min(displayedImage[i][c], min);
                max = std::max(displayedImage[i][c], max);
            });
            float d = max - min;
            if (d == 0) {
                displayedImage.iterateParallel([&](size_t i) {
                    displayedImage[i][c] = 0.f;
                });
            } else {
                displayedImage.iterateParallel([&](size_t i) {
                    displayedImage[i][c] = (displayedImage[i][c] - min) / d;
                });
            }
        }
    }

    if (img.isColor()) {
        displayedImage.iterateParallel([&](size_t i) {
            displayedImage[i] *= this->displayParameters.displayedColors;
        });
    } else {
        displayedImage.iterateParallel([&](size_t i) {
            displayedImage[i] = colorPalette(displayedImage[i].x(), this->displayParameters.colorRamp.getPath());
        });
    }
    return displayedImage;
}

QImage PlotImageData::computeDisplayedImage(const Vector3i& imgSize) const
{
    QImage emptyImg = QImage(imgSize.x(), imgSize.y(), QImage::Format_ARGB32); emptyImg.fill(Qt::white);
    if (this->image.empty()) return emptyImg;
    auto displayedImage = this->image.getColorImage(); //.resize(imgSize);
    if (displayedImage.empty()) return emptyImg;

    displayedImage = this->prepareImageForDisplay(this->image);
    unsigned char* data = new unsigned char[displayedImage.size() * 4];

    for (size_t i = 0; i < displayedImage.size(); ++i) {
        data[int(4 * i + 2)] = (unsigned char)(std::clamp(int(displayedImage[i].x() * 255), 0, 255));
        data[int(4 * i + 1)] = (unsigned char)(std::clamp(int(displayedImage[i].y() * 255), 0, 255));
        data[int(4 * i + 0)] = (unsigned char)(std::clamp(int(displayedImage[i].z() * 255), 0, 255));
        data[int(4 * i + 3)] = (unsigned char) 255;       // Alpha
    }

    QImage result = QImage(data, displayedImage.sizeX, displayedImage.sizeY, QImage::Format_ARGB32);
    if (imgSize.isValid()) result = result.scaled(imgSize.x(), imgSize.y());
    return result;
}

QImage PlotImageData::computeDisplayedImage(const GridV3 &overlay, const GridF &overlayAlpha) const
{
    return this->computeDisplayedImage({{"", overlay}}, {{"", overlayAlpha}}, {{"", true}}, {{"", 0.f}}, this->getImage().getDimensions());
}

QImage PlotImageData::computeDisplayedImage(const std::map<std::string, GridV3> &overlays,
                                            const std::map<std::string, GridF> &overlayAlphas,
                                            const std::map<std::string, bool>& displayedOverlays,
                                            const std::map<std::string, int>& overlayLayers,
                                            const Vector3i &imgSize) const
{
    Vector3i largestDimensions = imgSize;
    for (auto& [name, over] : overlays) {
        largestDimensions.x() = std::max(largestDimensions.x(), (int)over.sizeX);
        largestDimensions.y() = std::max(largestDimensions.y(), (int)over.sizeY);
    }
    if (!this->image.empty()) {
        largestDimensions.x() = std::max(this->image.getBwImage().getDimensions().x(), largestDimensions.x());
        largestDimensions.y() = std::max(this->image.getBwImage().getDimensions().y(), largestDimensions.y());
    }
    largestDimensions = Vector3i(largestDimensions.x(), largestDimensions.y(), 1);
    QImage img = this->computeDisplayedImage(largestDimensions);
    QPainter painter = QPainter(&img);

    auto sort = [=](std::map<std::string, int> M) -> std::vector<std::pair<std::string, int>> {

        // Declare vector of pairs
        std::vector<std::pair<std::string, int> > A;
        for (auto& it : M) {
            A.push_back(it);
        }
        std::sort(A.begin(), A.end(), [=](std::pair<std::string, int>& a, std::pair<std::string, int>& b){ return a.second < b.second; });
        return A;
    };

    auto overlayOrder = sort(overlayLayers);

    for (auto& [name, layerPriority] : overlayOrder) {
        if (displayedOverlays.count(name) && displayedOverlays.at(name) == false) continue;
        const auto& overlay = overlays.at(name); //.resize(imgSize, RESIZE_MODE::LINEAR);
        const auto& overlayAlpha = overlayAlphas.at(name); //.resize(imgSize, RESIZE_MODE::NEAREST);
        unsigned char* data = new unsigned char[overlay.size() * 4];

        for (size_t i = 0; i < overlay.size(); ++i) {
            data[int(4 * i + 2)] = (unsigned char)(std::clamp(int(overlay[i].x() * 255), 0, 255));
            data[int(4 * i + 1)] = (unsigned char)(std::clamp(int(overlay[i].y() * 255), 0, 255));
            data[int(4 * i + 0)] = (unsigned char)(std::clamp(int(overlay[i].z() * 255), 0, 255));
            data[int(4 * i + 3)] = (unsigned char) int((overlayAlpha.size() == overlay.size() ? overlayAlpha[i] : 1.f) * 255.f);       // Alpha
        }
        // std::cout << "Painting " << name << " at resolution " << overlay.sizeX << "x" << overlay.sizeY << " scaled to " << imgSize.x() << "x" << imgSize.y() << std::endl;
        painter.drawImage(0, 0, QImage(data, overlay.sizeX, overlay.sizeY, QImage::Format_ARGB32).scaled(largestDimensions.x(), largestDimensions.y()));
    }
    painter.end();
    if (!imgSize.isValid()) return img;
    return img.scaled(imgSize.x(), imgSize.y());
}

// void PlotImageData::setOnImageModified(const std::function<void ()> &callback)
// {
    // this->onImageModifiedCallbacks.push_back(callback);
// }
