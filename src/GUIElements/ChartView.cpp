#include "ChartView.h"


ChartView::ChartView(QWidget *parent) : QChartView(new Chart(nullptr), parent)
{
    // QObject::connect(this->chart(), &QChart::geometryChanged, this, &ChartView::updateLabelsPositions);
    // QObject::connect(this->chart(), &QChart::plotAreaChanged, this, &ChartView::updateLabelsPositions);
    // QObject::connect(this, &ChartView::updated, this, &ChartView::updateLabelsPositions);
}

ChartView& ChartView::setPlotModel(std::shared_ptr<PlotModel> dataModel, const std::string& title)
{
    this->_dataModel = std::move(dataModel);
    this->chart()->removeAllSeries();

    if (!title.empty())
        this->chart()->setTitle(QString::fromStdString(title));

    for (auto& labels : this->_dataModel->graphicLabels)
        for (auto& lab : labels)
            delete lab;
    this->_dataModel->graphicLabels.clear();

    int width = static_cast<int>(this->chart()->plotArea().width());
    int height = static_cast<int>(this->chart()->plotArea().height());
    int ViewW = static_cast<int>(this->width());
    int ViewH = static_cast<int>(this->height());
    QImage scaledImage = QImage(width, height, QImage::Format_ARGB32);
    scaledImage.fill(Qt::white);

    Vector3i renderResolution = Vector3i::invalid; // Vector3i(100, 100, 1); //Vector3i(clamp(width, 20, 400), clamp(height, 20, 400), 1);

    auto overlays = this->overlayColors;
    auto overlayAlphas = this->overlayAlpha;
    auto overlayDisplays = this->overlayDisplayed;
    auto overlayLayers = this->overlayLayer;
    if (!this->_dataModel->vectorData.field.empty()) {
        auto [overlay, alpha] = this->_dataModel->vectorData.getFieldImageAndAlpha(renderResolution, Vector3i(20, 20, 1));
        overlays["vector field"] = overlay;
        overlayAlphas["vector field"] = alpha;
        overlayDisplays["vector field"] = true;
        overlayLayers["vector field"] = 10000;
    }
    scaledImage = this->_dataModel->imageData.computeDisplayedImage(overlays, overlayAlphas, overlayDisplays, overlayLayers, renderResolution);
    scaledImage = scaledImage.scaled(QSize(width, height), Qt::IgnoreAspectRatio, Qt::TransformationMode::SmoothTransformation); // FastTransformation); // SmoothTransformation);

    QImage translated(ViewW, ViewH, QImage::Format_ARGB32);
    translated.fill(Qt::white);
    QPainter painter(&translated);
    QPointF TopLeft = this->chart()->plotArea().topLeft();
    painter.drawImage(TopLeft, scaledImage);

    //Display image in background
    //        this->chart()->setPlotAreaBackgroundBrush(scaledImage);
    this->chart()->setPlotAreaBackgroundBrush(translated);
    this->chart()->setPlotAreaBackgroundVisible(true);

    for (size_t i = 0; i < this->_dataModel->plot_data.size(); i++) {
        QLineSeries *series = new QLineSeries();
        if (this->_dataModel->plot_names.size() > 0 && this->_dataModel->plot_names.size() == this->_dataModel->plot_data.size())
            series->setName(QString::fromStdString(this->_dataModel->plot_names[i]));
        for (auto pos : this->_dataModel->plot_data[i])
            series->append(pos.x(), pos.y());
        //        if (PlotColorToQColor.count(this->_dataModel->plot_colors[i]))
        //            series->setColor(PlotColorToQColor.at(this->_dataModel->plot_colors[i]));
        series->setColor(this->_dataModel->plot_colors[i]);
        this->chart()->addSeries(series);
        if (series->name().isEmpty()) {
            this->chart()->legend()->markers(series)[0]->setVisible(false);
        }
    }
    for (size_t i = 0; i < this->_dataModel->scatter_data.size(); i++) {
        QScatterSeries *series = new QScatterSeries();
        if (this->_dataModel->scatter_names.size() > 0 && this->_dataModel->scatter_names.size() == this->_dataModel->scatter_data.size())
            series->setName(QString::fromStdString(this->_dataModel->scatter_names[i]));
        for (size_t j = 0; j < this->_dataModel->scatter_data[i].size(); j++) { //(auto pos : this->_dataModel->scatter_data[i]) {
            auto pos = this->_dataModel->scatter_data[i][j];
            series->append(pos.x(), pos.y());
        }
        this->chart()->addSeries(series);

        if (series->name().isEmpty()) {
            this->chart()->legend()->markers(series)[0]->setVisible(false);
        }
    }
    for (size_t iScatter = 0; iScatter < this->_dataModel->scatter_labels.size(); iScatter++) {
        this->_dataModel->graphicLabels.push_back(std::vector<QGraphicsTextItem*>());
        for (size_t iPoint = 0; iPoint < this->_dataModel->scatter_labels[iScatter].size(); iPoint++) {
            QGraphicsTextItem *itm = new QGraphicsTextItem(QString::fromStdString(this->_dataModel->scatter_labels[iScatter][iPoint]), this->chart());
            this->_dataModel->graphicLabels[iScatter].push_back(itm);
        }
    }

    this->chart()->createDefaultAxes();
    return *this;
}

ChartView& ChartView::updateLabelsPositions()
{
    if (this->_dataModel == nullptr) return *this;
    //    this->blockSignals(true);
    if (!this->_dataModel->selectedPlotData.empty() || !this->_dataModel->selectedScatterData.empty()) {
        QPointF qNewPoint = this->chart()->mapToValue(this->previousMousePos);
        Vector3 newPoint = Vector3(qNewPoint.x(), qNewPoint.y());
        for (auto& [iPlot, iPoint] : this->_dataModel->selectedPlotData)
            this->_dataModel->plot_data[iPlot][iPoint] = newPoint;
        for (auto& [iPlot, iPoint] : this->_dataModel->selectedScatterData)
            this->_dataModel->scatter_data[iPlot][iPoint] = newPoint;
    }

    for (size_t iScatter = 0; iScatter < this->_dataModel->scatter_labels.size(); iScatter++) {
        for (size_t iPoint = 0; iPoint < this->_dataModel->scatter_labels[iScatter].size(); iPoint++) {
            //                this->_dataModel->graphicLabels[iScatter][iPoint]->setPos(QPointF(this->_dataModel->scatter_data[iScatter][iPoint].first, this->_dataModel->scatter_data[iScatter][iPoint].second)); // this->chart()->mapToPosition(QPointF(this->scatter_data[iScatter][iPoint].first, this->scatter_data[iScatter][iPoint].second)));
            this->_dataModel->graphicLabels[iScatter][iPoint]->setPos(this->chart()->mapToPosition(QPointF(this->_dataModel->scatter_data[iScatter][iPoint].x(), this->_dataModel->scatter_data[iScatter][iPoint].y())));
        }
    }
    if (!this->_dataModel->selectedPlotData.empty() || !this->_dataModel->selectedScatterData.empty()) {
        // this->draw();
        std::cout << "Removed a call to ImageViewer::draw() here..." << std::endl;
    }
    //    this->blockSignals(false);
    return *this;
}

bool ChartView::selectData(const Vector3& pos)
{
    if (!pos.isValid()) return false;

    float minDist = 0.05f;
    this->_dataModel->selectedScatterData.clear();
    this->_dataModel->selectedPlotData.clear();

    if (pos.isValid()) {
        for (size_t i = 0; i < this->_dataModel->plot_data.size(); i++) {
            for (size_t j = 0; j < this->_dataModel->plot_data[i].size(); j++) {
                if ((this->_dataModel->plot_data[i][j] - pos).norm2() < minDist*minDist)
                    this->_dataModel->selectedPlotData.push_back({i, j});
            }
        }
        for (size_t i = 0; i < this->_dataModel->scatter_data.size(); i++) {
            for (size_t j = 0; j < this->_dataModel->scatter_data[i].size(); j++) {
                if ((this->_dataModel->scatter_data[i][j] - pos).norm2() < minDist*minDist)
                    this->_dataModel->selectedScatterData.push_back({i, j});
            }
        }
    }

    if (!this->_dataModel->selectedPlotData.empty() || !this->_dataModel->selectedScatterData.empty()) {
        this->lockView();
        return false;
    } else {
        this->unlockView();
        return true;
    }

    // if (!this->_dataModel->displayedImage.empty()) {
    // Q_EMIT clickedOnImage(pos * this->_dataModel->displayedImage.getDimensions(), this->_dataModel->displayedImage(pos * this->_dataModel->displayedImage.getDimensions()));
    // }
    // return *this;
}

Vector3 ChartView::getRelativeMousePositionInImage(const Vector3 &pos)
{
    // Get the coordinate in the plotted area...
    QPointF qMousePos(pos.x(), pos.y());
    QPointF mousePosInChart = this->chart()->mapFromParent(qMousePos);
    QRectF plotArea = this->chart()->plotArea();
    QPointF mousePosInPlot = mousePosInChart - plotArea.topLeft();
    QPointF qRelativeMousePos = mousePosInPlot;
    Vector3 mousePos(qRelativeMousePos.x() / float(plotArea.width()), qRelativeMousePos.y() / float(plotArea.height()));
    return mousePos;
}

const ChartView& ChartView::setOverlay(const GridV3 &image, std::string layerName, const GridF &alpha, int overlayLayer)
{
    this->overlayColors[layerName] = image;
    this->overlayAlpha[layerName] = alpha;
    this->overlayDisplayed[layerName] = true;
    this->overlayLayer[layerName] = overlayLayer;
    return *this;
}

bool ChartView::viewportEvent(QEvent *event)
{
    return QChartView::viewportEvent(event);
}

void ChartView::mousePressEvent(QMouseEvent *event)
{
    this->previousMousePos = event->pos();

    QPointF local = this->chart()->mapToValue(this->previousMousePos);
    Vector3 vecLocal = Vector3(local.x(), local.y());

    vecLocal = this->getRelativeMousePositionInImage(Vector3(event->pos().x(), event->pos().y()));
    this->selectData(vecLocal);
    Q_EMIT this->clickedOnValue(vecLocal, event->button() == Qt::LeftButton, event->button() == Qt::RightButton);

    QChartView::mousePressEvent(event);
}
void ChartView::mouseMoveEvent(QMouseEvent *event)
{
    if (event->buttons().testFlag(Qt::LeftButton)) {
        QPoint delta = event->pos() - this->previousMousePos;
        if (!this->locked)
            this->chart()->scroll(-delta.x(), delta.y());
        Q_EMIT this->updated();
    }

    auto mousePos = this->getRelativeMousePositionInImage(Vector3(event->pos().x(), event->pos().y()));
    auto prevPos = this->getRelativeMousePositionInImage(Vector3(previousMousePos.x(), previousMousePos.y()));


    this->previousMousePos = event->pos();
    Q_EMIT this->mouseMoved(mousePos, prevPos, event);
    return QChartView::mouseMoveEvent(event);
}
void ChartView::mouseReleaseEvent(QMouseEvent *event)
{
    this->selectData(Vector3::invalid);
    Q_EMIT this->clickedOnValue(Vector3::invalid, event->button() == Qt::LeftButton, event->button() == Qt::RightButton); // "unclick"
    return QChartView::mouseReleaseEvent(event);
}
void ChartView::keyPressEvent(QKeyEvent *event)
{
    switch (event->key()) {
    case Qt::Key_Plus:
        chart()->zoomIn();
        Q_EMIT this->updated();
        break;
    case Qt::Key_Minus:
        chart()->zoomOut();
        Q_EMIT this->updated();
        break;
    }
    Q_EMIT this->keyPressed(event);
    return QChartView::keyPressEvent(event);
}
void ChartView::wheelEvent(QWheelEvent *event)
{
    qreal factor = event->angleDelta().y() > 0? 0.5: 2.0;
    chart()->zoom(factor);
    Q_EMIT this->updated();
    event->accept();
    QChartView::wheelEvent(event);
}

Chart::Chart(QGraphicsItem* parent) : QChart(parent)
{}

bool Chart::sceneEvent(QEvent *event)
{
    return QChart::event(event);
}

bool Chart::gestureEvent(QGestureEvent *event)
{
    return true;
}


PlotModel::PlotModel()
{

}

PlotModel& PlotModel::addPlot(const std::vector<Vector3>& data, const std::string& name, const QColor& color)
{
    this->plot_data.push_back(data);
    this->plot_names.push_back(name);
    this->plot_colors.push_back(color);
    return *this;
}

PlotModel& PlotModel::addScatter(const std::vector<Vector3>& data, const std::string& name, const std::vector<std::string>& labels, std::vector<QColor> colors)
{
    if (colors.size() == 0) {
        colors = std::vector<QColor>({Qt::blue});
    }
    if (colors.size() == 1) {
        colors = std::vector<QColor>(data.size(), colors.front());
    }
    this->scatter_data.push_back(data);
    this->scatter_names.push_back(name);
    this->scatter_labels.push_back(labels);
    this->scatter_colors.push_back(colors);
    return *this;
}

PlotModel& PlotModel::addImage(const GridV3& image/*, bool clamped, bool normalized, bool absolute, const Vector3 &minColors, const Vector3 &maxColors*/)
{
    this->imageData.setImage(image);
    return *this;
}

PlotModel& PlotModel::addImage(const GridF &image)
{
    this->imageData.setImage(image);
    return *this;
}

PlotModel& PlotModel::addVectorField(const GridV3 &field)
{
    this->vectorData.setField(field);
    return *this;
}

PlotModel& PlotModel::reset()
{
    // this->backImage = nullptr;
    this->title = "";
    this->plot_data.clear();
    this->plot_names.clear();
    this->plot_colors.clear();
    this->scatter_data.clear();
    this->scatter_labels.clear();
    this->scatter_colors.clear();
    this->scatter_names.clear();
    this->graphicLabels.clear();

    this->selectedScatterData.clear();
    this->selectedPlotData.clear();

    this->imageData = PlotImageData();
    this->vectorData = PlotVectorData();
    return *this;
}

PlotImageData::PlotImageData() : PlotImageData(GridV3())
{

}

PlotImageData::PlotImageData(const GridV3 &img) : image(img)
{

}

PlotImageData::PlotImageData(const GridF &img) : image(img)
{

}

PlotImageData& PlotImageData::setImage(const GridV3 &img)
{
    this->image.setImage(img);
    return *this;
}

PlotImageData& PlotImageData::setImage(const GridF &img)
{
    this->image.setImage(img);
    return *this;
}

PlotImageData& PlotImageData::setNormalized(bool normalize)
{
    // this->normalized = normalize;
    this->displayParameters.normalized = normalize;
    return *this;
}

PlotImageData& PlotImageData::setColorRanges(const Vector3 &minRange, const Vector3 &maxRange)
{
    // this->colorRangeMin = minRange;
    // this->colorRangeMax = maxRange;
    this->displayParameters.colorRangeMin = minRange;
    this->displayParameters.colorRangeMax = maxRange;
    return *this;
}

PlotImageData& PlotImageData::setAbsolute(bool absolute)
{
    // this->absolute = absolute;
    this->displayParameters.absolute = absolute;
    return *this;
}

PlotImageData& PlotImageData::setClamped(bool clamp)
{
    // this->clamped = clamp;
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
            displayedImage[i] = colorPalette(displayedImage[i].x(), this->displayParameters.colorRamp.points);
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
        data[int(4 * i + 2)] = (unsigned char)(std::clamp(displayedImage[i].x(), 0.f, 1.f) * 255);
        data[int(4 * i + 1)] = (unsigned char)(std::clamp(displayedImage[i].y(), 0.f, 1.f) * 255);
        data[int(4 * i + 0)] = (unsigned char)(std::clamp(displayedImage[i].z(), 0.f, 1.f) * 255);
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
            data[int(4 * i + 2)] = (unsigned char)(std::clamp(overlay[i].x(), 0.f, 1.f) * 255);
            data[int(4 * i + 1)] = (unsigned char)(std::clamp(overlay[i].y(), 0.f, 1.f) * 255);
            data[int(4 * i + 0)] = (unsigned char)(std::clamp(overlay[i].z(), 0.f, 1.f) * 255);
            data[int(4 * i + 3)] = (unsigned char) int((overlayAlpha.size() == overlay.size() ? overlayAlpha[i] : 1.f) * 255.f);       // Alpha
        }
        // std::cout << "Painting " << name << " at resolution " << overlay.sizeX << "x" << overlay.sizeY << " scaled to " << imgSize.x() << "x" << imgSize.y() << std::endl;
        painter.drawImage(0, 0, QImage(data, overlay.sizeX, overlay.sizeY, QImage::Format_ARGB32).scaled(largestDimensions.x(), largestDimensions.y()));
    }
    painter.end();
    if (!imgSize.isValid()) return img;
    return img.scaled(imgSize.x(), imgSize.y());
}

PlotVectorData::PlotVectorData() : PlotVectorData(GridV3())
{}

PlotVectorData::PlotVectorData(const GridV3 &field) : field(field)
{}

PlotVectorData& PlotVectorData::setField(const GridV3 &field)
{
    this->field = field;
    return *this;
}

std::pair<GridV3, GridF> PlotVectorData::getFieldImageAndAlpha(const Vector3i &imgSize, const Vector3i& numberOfCells) const
{
    return PlotVectorData::createFieldImageAndAlpha(this->field, imgSize, numberOfCells, this->displayParameters);
}

std::pair<GridV3, GridF> PlotVectorData::createFieldImageAndAlpha(const GridV3 &field, Vector3i imgSize, const Vector3i &numberOfCells, DisplayedVectorFieldParameters displayParameters)
{
    const Vector3& backgroundColor = displayParameters.backgroundColor;
    if (!imgSize.isValid()) imgSize = field.getDimensions();
    GridV3 img(imgSize.x(), imgSize.y(), 1, backgroundColor);
    GridF alpha(img.getDimensions());
    if (field.empty()) return {img, alpha};

    GridV3 reduced = field.resize(numberOfCells);

    float minMag = std::numeric_limits<float>::max();
    float maxMag = std::numeric_limits<float>::lowest();
    reduced.iterate([&] (size_t i) {
        float mag = reduced[i].norm2();
        minMag = std::min(minMag, mag);
        maxMag = std::max(maxMag, mag);
    });
    minMag = std::sqrt(minMag);
    maxMag = std::sqrt(maxMag);
    // std::cout << minMag << " " << maxMag << std::endl;
    Vector3 reducedSize = reduced.getDimensions(); //imgSize / numberOfCells;
    Vector3 imageToReducedRatio = Vector3((float)imgSize.x() / (float)reducedSize.x(), (float)imgSize.y() / (float)reducedSize.y(), 1.f);
    Vector3 fieldToReducedRatio = Vector3((float)field.sizeX / (float)reducedSize.x(), (float)field.sizeY / (float)reducedSize.y(), 1.f);
    Vector3 fieldToImageRatio = Vector3((float)field.sizeX / (float)imgSize.x(), (float)field.sizeY / (float)imgSize.y(), 1.f);

    if (displayParameters.displayMode == DisplayedVectorFieldParameters::ARROWS) {
        reduced.iterateParallel([&] (const Vector3& _p) {
            Vector3 p = _p + Vector3(.5f, .5f);
            // AABBox cell((p - Vector3(.5f, .5f, 1)) * ratio, (p + Vector3(.5f, .5f, 1)) * ratio); // Added an depth (z) to avoid issue on the intersection computation
            Vector3 vec = reduced.interpolate(p);
            if (!vec.isValid()) return;
            float mag = vec.norm();
            if (mag < 1e-5) return;
            Vector3 dir = vec / mag;
            Vector3 color = displayParameters.colorRamp.getPoint(0.f);
            float relativeMag = 1.f;
            if (std::abs(minMag - maxMag) > 1e-5) {
                relativeMag = interpolation::linear(mag, minMag, maxMag);
                color = colorPalette(relativeMag, displayParameters.colorRamp.points);
            }
            bool valid = dir.xy().norm2() > 1e-5;
            if (!valid) return;

            Vector3 startLine = (p - dir * interpolation::inv_linear(relativeMag, .5f, 1.f)) * imageToReducedRatio;
            Vector3 endLine = (p + dir * interpolation::inv_linear(relativeMag, .5f, 1.f))  * imageToReducedRatio;
            float length = (endLine - startLine).norm();

            img = PlotVectorData::drawLine(img, color, startLine, endLine);
            alpha = PlotVectorData::drawLine(alpha, 1.f, startLine, endLine);

            img = PlotVectorData::drawLine(img, color, endLine, endLine - dir.rotated(deg2rad(20), Vector3(0, 0, 1)) * length * .3f);
            alpha = PlotVectorData::drawLine(alpha, 1.f, endLine, endLine - dir.rotated(deg2rad(20), Vector3(0, 0, 1)) * length * .3f);

            img = PlotVectorData::drawLine(img, color, endLine, endLine - dir.rotated(deg2rad(-20), Vector3(0, 0, 1)) * length * .3f);
            alpha = PlotVectorData::drawLine(alpha, 1.f, endLine, endLine - dir.rotated(deg2rad(-20), Vector3(0, 0, 1)) * length * .3f);
        });
    }
    else if (displayParameters.displayMode == DisplayedVectorFieldParameters::FLOWLINES) {
        int trailLength = 100;
        float stepLength = .1f;

        for (int x = 0; x < numberOfCells.x() - 1; x++) {
            for (int y = 0; y < numberOfCells.y() - 1; y++) {
                Vector3 p = Vector3(x + .5f, y + .5f) * fieldToReducedRatio;
                for (int i = 0; i < trailLength; i++) {
                    Vector3 dir = field.interpolate(p);
                    if (!dir.isValid()) break;
                    float mag = dir.norm();
                    if (mag < 1e-5) break;
                    auto color = colorPalette(interpolation::inv_linear(mag, minMag, maxMag), displayParameters.colorRamp.points);

                    dir = (dir.maxMagnitude(1.f)) * fieldToReducedRatio * stepLength;
                    Vector3 end = p + dir;

                    img = PlotVectorData::drawLine(img, color, p / fieldToImageRatio, end / fieldToImageRatio);
                    alpha = PlotVectorData::drawLine(alpha, 1.f, p / fieldToImageRatio, end / fieldToImageRatio);

                    p = end;
                    if (x == 0 && y == 0) std::cout << i << " -> " << p << "(" << mag << " / " << minMag << " / " << maxMag << ")" << std::endl;
                }
            }
        }
    }
    return {img, alpha};
}












