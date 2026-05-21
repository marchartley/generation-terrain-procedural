#include "ChartView.h"


ChartView::ChartView(QWidget *parent) : QChartView(new Chart(nullptr), parent)
{
}

ChartView& ChartView::setPlotModel(std::shared_ptr<PlotModel> dataModel, const std::string& title)
{
    this->_dataModel = std::move(dataModel);

    resetPlot();

    if (!title.empty())
        this->chart()->setTitle(QString::fromStdString(title));

    displayImages();

    this->chart()->setPlotAreaBackgroundVisible(true);

    displayPlotLines();
    displayScatterPoints();

    if (axisX == nullptr) {
        axisX = new QValueAxis(chart());
        axisY = new QValueAxis(chart());

        chart()->addAxis(axisX, Qt::AlignBottom);
        chart()->addAxis(axisY, Qt::AlignLeft);
    }

    Vector3 mini = Vector3::max();
    Vector3 maxi = Vector3::min();
    for (auto* s : chart()->series()) {
        s->attachAxis(axisX);
        s->attachAxis(axisY);

        if (auto series = qobject_cast<QXYSeries*>(s)) {
            for (const auto& p : series->points()) {
                mini = Vector3::min(mini, Vector3(p.x(), p.y()));
                maxi = Vector3::max(maxi, Vector3(p.x(), p.y()));
            }
        }
    }

    if (plottingLimits.dimensions().norm2() == 0) {
        plottingLimits = AABBox(mini, maxi);
        plottingLimits.expand({plottingLimits.min() - plottingLimits.dimensions() * .1f, plottingLimits.max() + plottingLimits.dimensions() * .1f});
    }

    axisX->setRange(plottingLimits.min().x(), plottingLimits.max().x());
    axisY->setRange(plottingLimits.min().y(), plottingLimits.max().y());

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
            plotLineData[iPlot][iPoint].pos = newPoint;
        for (auto& [iPlot, iPoint] : this->_dataModel->selectedScatterData)
            scatterData[iPlot][iPoint].pos = newPoint;
    }

    for (auto& series : scatterData) {
        for (auto& point : series) {
            auto graphicLabel = new QGraphicsTextItem(QString::fromStdString(point.label), this->chart());
            graphicLabel->setPos(point.pos.x(), point.pos.y());
        }
    }
    if (!this->_dataModel->selectedPlotData.empty() || !this->_dataModel->selectedScatterData.empty()) {
        std::cout << "Removed a call to ImageViewer::draw() here..." << std::endl;
    }
    return *this;
}

ChartView &ChartView::setPlottingLimits(const Vector3 &mini, const Vector3 &maxi)
{
    if (maxi.isValid()) {
        this->plottingLimits = AABBox(mini, maxi);
    } else {
        this->plottingLimits = AABBox(Vector3::origin, mini);
    }
    return *this;
}

bool ChartView::selectData(const Vector3& pos)
{
    if (!pos.isValid()) return false;

    auto& scatterData = _dataModel->scatterData;
    auto& plotLineData = _dataModel->plotLineData;

    float minDist = (this->axisX->max() - this->axisX->min()) * 0.05f;
    // float minDist = 0.05f;
    this->_dataModel->selectedScatterData.clear();
    this->_dataModel->selectedPlotData.clear();

    if (pos.isValid()) {
        for (size_t i = 0; i < plotLineData.size(); i++) {
            for (size_t j = 0; j < plotLineData[i].size(); j++) {
                if ((plotLineData[i][j].pos - pos).norm2() < minDist*minDist)
                    this->_dataModel->selectedPlotData.push_back({i, j});
            }
        }
        for (size_t i = 0; i < scatterData.size(); i++) {
            for (size_t j = 0; j < scatterData[i].size(); j++) {
                if ((scatterData[i][j].pos - pos).norm2() < minDist*minDist)
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

void ChartView::resetPlot()
{
    auto& imageData = _dataModel->imageData;
    auto& vectorData = _dataModel->vectorData;
    auto& scatterData = _dataModel->scatterData;
    auto& plotLineData = _dataModel->plotLineData;

    this->chart()->removeAllSeries();

    // this->chart()->createDefaultAxes();
}

void ChartView::displayImages()
{
    auto& vectorData = _dataModel->vectorData;
    auto& imageData = _dataModel->imageData;

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
    else {
        this->chart()->setPlotAreaBackgroundBrush(QBrush());
    }
}

void ChartView::displayPlotLines()
{
    auto& plotLineData = _dataModel->plotLineData;
    for (auto& series : plotLineData) {
        QLineSeries *q_series = new QLineSeries();
        for (const auto& point : series)
            q_series->append(point.pos.x(), point.pos.y());
        QPen pen = q_series->pen();
        QColor col = QColor(int(series.color.r() * 255), int(series.color.g() * 255), int(series.color.b() * 255));
        pen.setColor(col);
        q_series->setPen(pen);
        q_series->setColor(col);
        this->chart()->addSeries(q_series);
        if (series.name != "") {
            q_series->setName(QString::fromStdString(series.name));
        } else {
            this->chart()->legend()->markers(q_series)[0]->setVisible(false);
        }
    }
}

void ChartView::displayScatterPoints()
{
    auto& scatterData = _dataModel->scatterData;
    for (auto& series : scatterData) {
        QScatterSeries *q_series = new QScatterSeries();
        for (auto& point : series) {
            q_series->append(point.pos.x(), point.pos.y());
        }
        QPen pen = q_series->pen();
        QColor col = QColor(int(series.color.r() * 255), int(series.color.g() * 255), int(series.color.b() * 255));
        pen.setColor(col);
        q_series->setPen(pen);
        q_series->setColor(col);
        this->chart()->addSeries(q_series);
        if (series.name != "") {
            q_series->setName(QString::fromStdString(series.name));
        } else {
            this->chart()->legend()->markers(q_series)[0]->setVisible(false);
        }
    }
}

ChartView& ChartView::useCurrentPlottingLimits()
{
    this->plottingLimits.mini = Vector3(axisX->min(), axisY->min());
    this->plottingLimits.maxi = Vector3(axisX->max(), axisY->max());
    return *this;
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
    Q_EMIT this->mousePressed(event);
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
    Q_EMIT this->mouseReleased(event);
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
    float zoomFactor = 1.1f;
    qreal factor = event->angleDelta().y() > 0? 1.0/zoomFactor : zoomFactor;
    chart()->zoom(factor);
    Q_EMIT this->updated();
    event->accept();
    QChartView::wheelEvent(event);
}

Chart::Chart(QGraphicsItem* parent) : QChart(parent)
{}

bool Chart::sceneEvent(QEvent *event)
{
    return QChart::sceneEvent(event);
}

bool Chart::gestureEvent(QGestureEvent *event)
{
    return true;
}


PlotModel::PlotModel()
{

}

PlotModel& PlotModel::addPlot(const std::vector<Vector3>& data, const std::string& name, const Vector3& color)
{
    this->plotLineData.add(data, name, color);
    return *this;
}

PlotModel& PlotModel::addScatter(const std::vector<Vector3>& data, const std::string& name, const std::vector<std::string>& labels, const Vector3& color)
{
    // if (colors.size() == 0) {
    //     colors = std::vector<Vector3>({Vector3::blue});
    // }
    // if (colors.size() == 1) {
    //     colors = std::vector<Vector3>(data.size(), colors.front());
    // }
    this->scatterData.add(data, labels, color, name);
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

    // this->selectedScatterData.clear();
    // this->selectedPlotData.clear();

    this->imageData = PlotImageData();
    this->vectorData = PlotVectorData();
    return *this;
}


