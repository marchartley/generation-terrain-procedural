#include "ChartView.h"


ChartView::ChartView(QWidget *parent) : QChartView(new Chart(nullptr), parent)
{}

ChartView& ChartView::setPlotModel(std::shared_ptr<PlotModel> dataModel, const std::string& title)
{
    this->_dataModel = std::move(dataModel);

    auto& imageData = _dataModel->imageData;
    auto& vectorData = _dataModel->vectorData;
    auto& scatterData = _dataModel->scatterData;
    auto& plotLineData = _dataModel->plotLineData;

    this->chart()->removeAllSeries();

    if (!title.empty())
        this->chart()->setTitle(QString::fromStdString(title));

    for (auto& labels : scatterData.graphicLabels)
        for (auto& lab : labels)
            delete lab;
    scatterData.graphicLabels.clear();

    if (!imageData.image.empty() || !vectorData.field.empty() || !this->overlayColors.empty()) {
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
        if (!vectorData.field.empty()) {
            auto [overlay, alpha] = vectorData.getFieldImageAndAlpha(renderResolution, Vector3i(20, 20, 1));
            overlays["vector field"] = overlay;
            overlayAlphas["vector field"] = alpha;
            overlayDisplays["vector field"] = true;
            overlayLayers["vector field"] = 10000;
        }
        scaledImage = imageData.computeDisplayedImage(overlays, overlayAlphas, overlayDisplays, overlayLayers, renderResolution);
        scaledImage = scaledImage.scaled(QSize(width, height), Qt::IgnoreAspectRatio, Qt::TransformationMode::SmoothTransformation); // FastTransformation); // SmoothTransformation);

        QImage translated(ViewW, ViewH, QImage::Format_ARGB32);
        translated.fill(Qt::white);
        QPainter painter(&translated);
        QPointF TopLeft = this->chart()->plotArea().topLeft();
        painter.drawImage(TopLeft, scaledImage);

        this->chart()->setPlotAreaBackgroundBrush(translated);
    }
    this->chart()->setPlotAreaBackgroundVisible(true);

    for (size_t i = 0; i < plotLineData.plot_data.size(); i++) {
        QLineSeries *series = new QLineSeries();
        if (plotLineData.plot_names.size() > 0 && plotLineData.plot_names.size() == plotLineData.plot_data.size())
            series->setName(QString::fromStdString(plotLineData.plot_names[i]));
        for (auto pos : plotLineData.plot_data[i])
            series->append(pos.x(), pos.y());
        series->setColor(plotLineData.plot_colors[i]);
        this->chart()->addSeries(series);
        if (series->name().isEmpty()) {
            this->chart()->legend()->markers(series)[0]->setVisible(false);
        }
    }
    for (size_t i = 0; i < scatterData.scatter_data.size(); i++) {
        QScatterSeries *series = new QScatterSeries();
        if (scatterData.scatter_names.size() > 0 && scatterData.scatter_names.size() == scatterData.scatter_data.size())
            series->setName(QString::fromStdString(scatterData.scatter_names[i]));
        for (size_t j = 0; j < scatterData.scatter_data[i].size(); j++) {
            auto pos = scatterData.scatter_data[i][j];
            series->append(pos.x(), pos.y());
        }
        this->chart()->addSeries(series);

        if (series->name().isEmpty()) {
            this->chart()->legend()->markers(series)[0]->setVisible(false);
        }
    }
    for (size_t iScatter = 0; iScatter < scatterData.scatter_labels.size(); iScatter++) {
        scatterData.graphicLabels.push_back(std::vector<QGraphicsTextItem*>());
        for (size_t iPoint = 0; iPoint < scatterData.scatter_labels[iScatter].size(); iPoint++) {
            QGraphicsTextItem *itm = new QGraphicsTextItem(QString::fromStdString(scatterData.scatter_labels[iScatter][iPoint]), this->chart());
            scatterData.graphicLabels[iScatter].push_back(itm);
        }
    }

    this->chart()->createDefaultAxes();
    return *this;
}

ChartView& ChartView::updateLabelsPositions()
{
    if (this->_dataModel == nullptr) return *this;

    auto& scatterData = _dataModel->scatterData;
    auto& plotLineData = _dataModel->plotLineData;

    if (!this->_dataModel->selectedPlotData.empty() || !this->_dataModel->selectedScatterData.empty()) {
        QPointF qNewPoint = this->chart()->mapToValue(this->previousMousePos);
        Vector3 newPoint = Vector3(qNewPoint.x(), qNewPoint.y());
        for (auto& [iPlot, iPoint] : this->_dataModel->selectedPlotData)
            plotLineData.plot_data[iPlot][iPoint] = newPoint;
        for (auto& [iPlot, iPoint] : this->_dataModel->selectedScatterData)
            scatterData.scatter_data[iPlot][iPoint] = newPoint;
    }

    for (size_t iScatter = 0; iScatter < scatterData.scatter_labels.size(); iScatter++) {
        for (size_t iPoint = 0; iPoint < scatterData.scatter_labels[iScatter].size(); iPoint++) {
            scatterData.graphicLabels[iScatter][iPoint]->setPos(this->chart()->mapToPosition(QPointF(scatterData.scatter_data[iScatter][iPoint].x(), scatterData.scatter_data[iScatter][iPoint].y())));
        }
    }
    if (!this->_dataModel->selectedPlotData.empty() || !this->_dataModel->selectedScatterData.empty()) {
        std::cout << "Removed a call to ImageViewer::draw() here..." << std::endl;
    }
    return *this;
}

bool ChartView::selectData(const Vector3& pos)
{
    if (!pos.isValid()) return false;

    auto& scatterData = _dataModel->scatterData;
    auto& plotLineData = _dataModel->plotLineData;

    float minDist = 0.05f;
    this->_dataModel->selectedScatterData.clear();
    this->_dataModel->selectedPlotData.clear();

    if (pos.isValid()) {
        for (size_t i = 0; i < plotLineData.plot_data.size(); i++) {
            for (size_t j = 0; j < plotLineData.plot_data[i].size(); j++) {
                if ((plotLineData.plot_data[i][j] - pos).norm2() < minDist*minDist)
                    this->_dataModel->selectedPlotData.push_back({i, j});
            }
        }
        for (size_t i = 0; i < scatterData.scatter_data.size(); i++) {
            for (size_t j = 0; j < scatterData.scatter_data[i].size(); j++) {
                if ((scatterData.scatter_data[i][j] - pos).norm2() < minDist*minDist)
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
    this->plotLineData.add(data, name, color);
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
    this->scatterData.add(data, labels, colors, name);
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
    plotLineData.reset();
    scatterData.reset();

    this->selectedScatterData.clear();
    this->selectedPlotData.clear();

    this->imageData = PlotImageData();
    this->vectorData = PlotVectorData();
    return *this;
}


