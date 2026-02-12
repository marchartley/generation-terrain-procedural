#include "AbstractPlotter.h"


#include <iostream>


ChartView::ChartView(QWidget *parent) : QChartView(nullptr, parent)
{}

ChartView::ChartView(QChart *chart, QWidget *parent) : QChartView(chart, parent)
{

}

ChartView::ChartView(Chart *chart, QWidget *parent) : QChartView((QChart*)chart, parent)
{
    this->_chart = chart;

    QObject::connect(this->chart(), &QChart::geometryChanged, this, &ChartView::updateLabelsPositions);
    QObject::connect(this->chart(), &QChart::plotAreaChanged, this, &ChartView::updateLabelsPositions);
    QObject::connect(this, &ChartView::updated, this, &ChartView::updateLabelsPositions);
}

ChartView* ChartView::setPlotModel(PlotModel *dataModel, std::string title)
{
    this->_dataModel = dataModel;
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
    // if (!this->_dataModel->getImage().empty()) {
        auto overlays = this->overlayColors;
        auto overlayAlphas = this->overlayAlpha;
        auto overlayDisplays = this->overlayDisplayed;
        if (!this->_dataModel->vectorData.field.empty()) {
            auto [overlay, alpha] = this->_dataModel->vectorData.getFieldImageAndAlpha(renderResolution, Vector3i(20, 20, 1));
            overlays["vector field"] = overlay;
            overlayAlphas["vector field"] = alpha;
            overlayDisplays["vector field"] = true;
        }
        //scale the image to fit plot area
        // if (!overlays.empty())
            scaledImage = this->_dataModel->imageData.computeDisplayedImage(overlays, overlayAlphas, overlayDisplays, renderResolution);
        // else
        // this->_dataModel->imageData.image.setImage(this->_dataModel->vectorData.getFieldImageAndAlpha(renderResolution, Vector3i(20, 20, 1)).first);
            // scaledImage = this->_dataModel->imageData.computeDisplayedImage(renderResolution);
    // }
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
    return this;
}

ChartView* ChartView::updateLabelsPositions()
{
    if (this->_dataModel == nullptr) return this;
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
    return this;
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
    // return this;
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

ChartView *ChartView::setOverlay(const GridV3 &image, std::string layerName, const GridF &alpha)
{
    this->overlayColors[layerName] = image;
    this->overlayAlpha[layerName] = alpha;
    this->overlayDisplayed[layerName] = true;
    return this;
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

Chart::Chart() : QChart()
{}

bool Chart::sceneEvent(QEvent *event)
{
    //    if (event->type() == QEvent::Gesture) {
    //        return this->gestureEvent(static_cast<QGestureEvent *>(event));
    //    }
    return QChart::event(event);
}

bool Chart::gestureEvent([[maybe_unused]] QGestureEvent *event)
{
    //    if (QGesture *gesture = event->gesture(Qt::PanGesture)) {
    //        QPanGesture *pan = static_cast<QPanGesture *>(gesture);
    //        QChart::scroll(-(pan->delta().x()), pan->delta().y());
    //    }

    //    if (QGesture *gesture = event->gesture(Qt::PinchGesture)) {
    //        QPinchGesture *pinch = static_cast<QPinchGesture *>(gesture);
    //        if (pinch->changeFlags() & QPinchGesture::ScaleFactorChanged)
    //            QChart::zoom(pinch->scaleFactor());
    //    }

    return true;
}

//void Chart::wheelEvent(QWheelEvent *event)
//{
//    QPoint scroll = event->angleDelta();
//    QChart::scroll(scroll.x(), scroll.y());
//    return QChart::wheelEvent(event);
//}




















std::map<std::string, AbstractPlotter*> AbstractPlotter::instances = std::map<std::string, AbstractPlotter*>();
std::string AbstractPlotter::defaultName = "default";

AbstractPlotter::AbstractPlotter(const std::string& name, const std::string& title, QWidget *parent) : AbstractPlotter(name, new ChartView(new Chart()), title, parent)
{
}

AbstractPlotter::AbstractPlotter(const std::string& name, ChartView *chartView, const std::string &title, QWidget *parent) : QDialog(parent), chartView(chartView), name(name)
{
    if (this->chartView == nullptr)
        this->chartView = new ChartView(new Chart());

    this->dataModel = new PlotModel();

    auto layout = new InterfaceUI(new QVBoxLayout());
    auto mainLayout = new InterfaceUI(new QHBoxLayout());

    //    auto right = new QVBoxLayout();
    mainInterface = new InterfaceUI(new QVBoxLayout(), "Main");
    mainInterface->add(new UIElement(this->chartView));
    toolsInterface = new InterfaceUI(new QVBoxLayout(), "Tools");
    viewOptionsInterface = new InterfaceUI(new QVBoxLayout(), "View options");
    saveCopyInterface = new InterfaceUI(new QVBoxLayout(), "Save/Copy");
    infosInterface = new InterfaceUI(new QHBoxLayout(), "Infos");

    this->chartView->setRenderHint(QPainter::Antialiasing);
    this->chartView->chart()->legend()->setMarkerShape(QLegend::MarkerShapeFromSeries);
    viewAndCopyInterface = new InterfaceUI(new QVBoxLayout());
    viewAndCopyInterface->add(std::vector<UIElement*>({viewOptionsInterface, saveCopyInterface}));
    //    this->chartView->setMaximumSize(10000, 10000);
    //    this->chartView->chart()->setMaximumSize(10000, 10000);
    //    this->chartView->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    //    this->chartView->chart()->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    //    this->mouseInfoLabel->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Maximum);

    //    this->setWindowModality(Qt::WindowModality::NonModal);
    //    this->setModal(false);

    layout->add(std::vector<UIElement*>({mainLayout, infosInterface}));

    mainLayout->add({toolsInterface, mainInterface, viewAndCopyInterface});

    this->viewAndCopyInterface->get()->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
    this->toolsInterface->get()->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
    this->mainInterface->get()->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
    qobject_cast<QHBoxLayout*>(mainLayout->get()->layout())->setStretchFactor(this->viewAndCopyInterface->get(), 1);
    qobject_cast<QHBoxLayout*>(mainLayout->get()->layout())->setStretchFactor(this->mainInterface->get(), 3);
    qobject_cast<QHBoxLayout*>(mainLayout->get()->layout())->setStretchFactor(this->toolsInterface->get(), 1);

    this->mouseInfoLabel = new QLabel("");
    infosInterface->add(new UIElement(this->mouseInfoLabel));

    auto saveButton = new ButtonElement("Save");
    auto copyToClipboardButton = new ButtonElement("Copy");
    saveButton->setOnClick([&]() {
        QString q_filename = QFileDialog::getSaveFileName(this, QString("Save plot"));
        if (!q_filename.isEmpty())
            saveFig(q_filename.toStdString());
    });
    copyToClipboardButton->setOnClick([&]() { this->copyToClipboard(); });

    saveCopyInterface->add(saveButton);
    saveCopyInterface->add(copyToClipboardButton);

    this->setLayout(layout->get()->layout());

    this->resize(1100, 600);
    this->updateGeometry();

    this->setWindowTitle(QString::fromStdString(toCapitalize(title)));

    // QObject::connect(this->chartView, &ChartView::clickedOnValue, this, &AbstractPlotter::selectData);
    QObject::connect(this->chartView, &ChartView::mouseMoved, this, [&](const Vector3& pos, const Vector3& prevPos, QMouseEvent* e){
        this->displayInfoUnderMouse(pos);
        if (this->hasImage()) {
            Q_EMIT this->movedOnImage(pos * dataModel->getImage().getDimensions(), prevPos * dataModel->getImage().getDimensions(), e);
        }
        else if (this->hasVectorField()) {
            Q_EMIT this->movedOnImage(pos * dataModel->vectorData.field.getDimensions(), prevPos * dataModel->vectorData.field.getDimensions(), e);
        }
    });
    QObject::connect(this->chartView->chart(), &QChart::geometryChanged, this, &AbstractPlotter::draw);
    QObject::connect(this->chartView, &ChartView::clickedOnValue, this, [&](const Vector3& pos, bool leftClick, bool rightClick) {
        if (this->hasImage()) {
            Q_EMIT this->clickedOnImage(pos, this->dataModel->getImage().at(pos * this->dataModel->getImage().getDimensions()), leftClick, rightClick);
        } else if (this->hasVectorField()) {
            Q_EMIT this->clickedOnImage(pos, Vector3::invalid, leftClick, rightClick);
        }
    });
}

AbstractPlotter* AbstractPlotter::addPlot(std::vector<float> data, std::string name, QColor color)
{
    std::vector<Vector3> _data;
    for (unsigned int i = 0; i < data.size(); i++) {
        _data.push_back(Vector3(i, data[i]));
    }
    return this->addPlot(_data, name, color);
}

AbstractPlotter* AbstractPlotter::addPlot(std::vector<Vector3> data, std::string name, QColor color)
{
    this->dataModel->addPlot(data, name, color);
    return this;
}

AbstractPlotter *AbstractPlotter::addPlot(const BSpline &data, std::string name, QColor color)
{
    return this->addPlot(data.points, name, color);
}

AbstractPlotter* AbstractPlotter::addScatter(std::vector<float> data, std::string name, std::vector<std::string> labels, std::vector<QColor> colors)
{
    std::vector<Vector3> _data;
    for (unsigned int i = 0; i < data.size(); i++) {
        _data.push_back(Vector3(i, data[i]));
    }
    return this->addScatter(_data, name, labels, colors);
}

AbstractPlotter* AbstractPlotter::addScatter(std::vector<Vector3> data, std::string name, std::vector<std::string> labels, std::vector<QColor> colors)
{
    this->dataModel->addScatter(data, name, labels, colors);
    return this;
}

AbstractPlotter* AbstractPlotter::addImage(const GridV3& image)
{
    this->dataModel->addImage(image);
    return this;
}

AbstractPlotter* AbstractPlotter::addImage(const GridF &image)
{
    this->dataModel->addImage(image);
    return this;
}

AbstractPlotter* AbstractPlotter::addImage(const Matrix3<double> &image)
{
    return this->addImage((GridF)image);
}

AbstractPlotter* AbstractPlotter::addImage(const GridI &image)
{
    return this->addImage((GridF)image);
}

AbstractPlotter *AbstractPlotter::addVectorField(const GridV3 &field)
{
    this->dataModel->addVectorField(field);
    return this;
}

AbstractPlotter* AbstractPlotter::setOverlay(const GridV3 &colors, const GridF &alpha, const std::string &overlayName)
{
    this->chartView->setOverlay(colors, overlayName, alpha);
    return this;
}

AbstractPlotter *AbstractPlotter::showOverlay(const std::string &overlayName)
{
    this->chartView->overlayDisplayed[overlayName] = true;
    return this;
}
AbstractPlotter *AbstractPlotter::hideOverlay(const std::string &overlayName)
{
    this->chartView->overlayDisplayed[overlayName] = false;
    return this;
}

GridV3 AbstractPlotter::computeStreamLinesRendering(const GridV3 &field, Vector3 imgSize) const
{
    if (!imgSize.isValid())
        imgSize = field.getDimensions();
    imgSize.z() = 1;
    GridV3 img(imgSize);
    GridF particlesPositions(20, 20, 1);
    Vector3 ratio = imgSize / field.getDimensions();

    int linesLength = 30;

    // Vector3 color(0, 1, 0);

    std::vector<Vector3> particles(particlesPositions.size());
    for (int i = 0; i < particles.size(); i++) {
        auto& particle = particles[i];
        particle = ((particlesPositions.getCoordAsVector3(i) + Vector3(.5, .5, 0)) / particlesPositions.getDimensions()) * field.getDimensions();
        Vector3 dir;

        BSpline spline;
        for (int t = 0; t < linesLength; t++) {
            dir = field.interpolate(particle).normalized();
            spline.points.push_back(particle);
            particle += dir;
        }
        auto path = spline.getPath(2.f * linesLength * ratio.maxComp());
        for (int t = 0; t < path.size(); t++)
            img(path[t] * ratio) = colorPalette(float(t) / float(path.size() - 1), Vector3(0, 1, 0), Vector3(.5, 1, 0));
    }
    return img;
}

AbstractPlotter* AbstractPlotter::addStreamLines(const GridV3 &field, Vector3 imgSize, float opacity)
{
    GridV3 img = computeStreamLinesRendering(field, imgSize);
    if (this->hasImage()) {
        img = this->dataModel->getImage().resize(img.getDimensions()) + img * (opacity);
    }
    return this->addImage(img);
}

AbstractPlotter* AbstractPlotter::draw()
{
    this->chartView->setPlotModel(this->dataModel);
    return this;
}

AbstractPlotter* AbstractPlotter::show()
{
    this->updateUI();
    this->draw();
    QDialog::show();
    return this;
}

AbstractPlotter* AbstractPlotter::updateUI(bool forceUpdate)
{
    if (forceUpdate || !this->isVisible()) {
        blockSignals(true);

        this->updateToolsInterface();
        this->updateViewOptionsInterface();
        this->updateSaveCopyInterface();
        this->updateInfosInterface();
        blockSignals(false);
    }
    return this;
}

AbstractPlotter *AbstractPlotter::updateToolsInterface()
{
    this->toolsInterface->update();
    return this;
}

AbstractPlotter *AbstractPlotter::updateViewOptionsInterface()
{
    this->viewOptionsInterface->update();
    return this;
}

AbstractPlotter *AbstractPlotter::updateSaveCopyInterface()
{
    this->saveCopyInterface->update();
    return this;
}

AbstractPlotter *AbstractPlotter::updateInfosInterface()
{
    this->infosInterface->update();
    return this;
}

int AbstractPlotter::exec()
{
    this->show();

    return QDialog::exec();

}

AbstractPlotter* AbstractPlotter::saveFig(std::string filename)
{
    QPixmap p = this->chartView->grab();
    if (this->hasImage())
        p = QPixmap::fromImage(this->dataModel->imageData.computeDisplayedImage());
    p.save(QString::fromStdString(filename), "PNG");
    std::cout << "Image " << filename << " saved." << std::endl;
    return this;
}

AbstractPlotter* AbstractPlotter::copyToClipboard()
{
    QPixmap p = this->chartView->grab();
    QApplication::clipboard()->setPixmap(p, QClipboard::Clipboard);
    return this;
}

void AbstractPlotter::resizeEvent(QResizeEvent *event)
{
    QDialog::resizeEvent(event);
    this->draw();
}

void AbstractPlotter::showEvent(QShowEvent *event)
{
    QDialog::showEvent(event);
    this->draw();
}

QTimer *AbstractPlotter::animate(std::function<void ()> callback, int interval_ms)
{
    QTimer* t = new QTimer(this);
    t->setInterval(interval_ms);
    t->setSingleShot(true);
    QObject::connect(t, &QTimer::timeout, [t, callback, interval_ms]() {
        callback();
        t->start(interval_ms);
    });
    t->start();
    return t;
}

AbstractPlotter* AbstractPlotter::reset()
{
    this->dataModel->reset();
    this->chartView->setPlotModel(this->dataModel);
    return this;
}

AbstractPlotter* AbstractPlotter::displayInfoUnderMouse(const Vector3 &relativeMousePos)
{
    if (!this->hasImage() || relativeMousePos.minComp() < 0.f || relativeMousePos.maxComp() > 1.f)
        return this;
    std::ostringstream oss;
    Vector3 size = this->dataModel->getImage().getDimensions();
    Vector3 position = relativeMousePos * size;
    Vector3 value = this->dataModel->getImage()(position);
    oss << "Mouse pos: " << int(position.x()) << ", " << int(position.y()) << " -- Value : (" << value.x() << ", " << value.y() << ", " << value.z() << ") ";
    this->mouseInfoLabel->setText(QString::fromStdString(oss.str()));
    return this;
}




TextItem::TextItem(QString text, QPoint pos, QChart *chart, QAbstractSeries *series)
    : QGraphicsItem(chart), _chart(chart), _series(series), _anchor(pos) {
    setText(text);
}

QRectF TextItem::boundingRect() const {
    QPointF anchor = mapFromParent(_chart->mapToPosition(_anchor, _series));
    QRectF rect;
    rect.setLeft(qMin(_textRect.left(), anchor.x()));
    rect.setRight(qMax(_textRect.right(), anchor.x()));
    rect.setTop(qMin(_textRect.top(), anchor.y()));
    rect.setBottom(qMax(_textRect.bottom(), anchor.y()));
    return rect;
}

void TextItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) {
    Q_UNUSED(option)
    Q_UNUSED(widget)
    QPointF anchor = mapFromParent(_chart->mapToPosition(_anchor, _series));
    painter->drawText(anchor, _text);
}

void TextItem::setText(const QString &text) {
    _text = text;
    QFontMetrics metrics((QFont()));
    _textRect = metrics.boundingRect(QRect(0, 0, 150, 150),
                                     Qt::AlignLeft, _text);
    _textRect.translate(5, 5);
    prepareGeometryChange();
}

void TextItem::setAnchor(QPointF point) { _anchor = point; }



PlotModel::PlotModel()
{

}

PlotModel *PlotModel::addPlot(const std::vector<Vector3>& data, const std::string& name, const QColor& color)
{
    this->plot_data.push_back(data);
    this->plot_names.push_back(name);
    this->plot_colors.push_back(color);
    return this;
}

PlotModel *PlotModel::addScatter(const std::vector<Vector3>& data, const std::string& name, const std::vector<std::string>& labels, std::vector<QColor> colors)
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
    return this;
}

PlotModel *PlotModel::addImage(const GridV3& image/*, bool clamped, bool normalized, bool absolute, const Vector3 &minColors, const Vector3 &maxColors*/)
{
    this->imageData.setImage(image);
    return this;
}

PlotModel *PlotModel::addImage(const GridF &image)
{
    this->imageData.setImage(image);
    return this;
}

PlotModel *PlotModel::addVectorField(const GridV3 &field)
{
    this->vectorData.setField(field);
    return this;
}

PlotModel* PlotModel::reset()
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
    return this;
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

PlotImageData *PlotImageData::setImage(const GridV3 &img)
{
    this->image.setImage(img);
    return this;
}

PlotImageData *PlotImageData::setImage(const GridF &img)
{
    this->image.setImage(img);
    return this;
}

PlotImageData *PlotImageData::setNormalized(bool normalize)
{
    // this->normalized = normalize;
    this->displayParameters.normalized = normalize;
    return this;
}

PlotImageData *PlotImageData::setColorRanges(const Vector3 &minRange, const Vector3 &maxRange)
{
    // this->colorRangeMin = minRange;
    // this->colorRangeMax = maxRange;
    this->displayParameters.colorRangeMin = minRange;
    this->displayParameters.colorRangeMax = maxRange;
    return this;
}

PlotImageData *PlotImageData::setAbsolute(bool absolute)
{
    // this->absolute = absolute;
    this->displayParameters.absolute = absolute;
    return this;
}

PlotImageData *PlotImageData::setClamped(bool clamp)
{
    // this->clamped = clamp;
    this->displayParameters.clamped = clamp;
    return this;
}

QImage PlotImageData::computeDisplayedImage(const Vector3i& imgSize) const
{
    QImage emptyImg = QImage(imgSize.x(), imgSize.y(), QImage::Format_ARGB32); emptyImg.fill(Qt::white);
    if (this->image.empty()) return emptyImg;
    auto displayedImage = this->image.getColorImage(); //.resize(imgSize);
    if (displayedImage.empty()) return emptyImg;

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

    if (this->image.isColor()) {
        displayedImage.iterateParallel([&](size_t i) {
            displayedImage[i] *= this->displayParameters.displayedColors;
        });
    } else {
        displayedImage.iterateParallel([&](size_t i) {
            displayedImage[i] = colorPalette(displayedImage[i].x(), this->displayParameters.colorRamp.points);
        });
    }

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
    /*QImage img = this->computeDisplayedImage();
    QPainter painter = QPainter(&img);

    unsigned char* data = new unsigned char[overlay.size() * 4];

    for (size_t i = 0; i < overlay.size(); ++i) {
        data[int(4 * i + 2)] = (unsigned char)(std::clamp(overlay[i].x(), 0.f, 1.f) * 255);
        data[int(4 * i + 1)] = (unsigned char)(std::clamp(overlay[i].y(), 0.f, 1.f) * 255);
        data[int(4 * i + 0)] = (unsigned char)(std::clamp(overlay[i].z(), 0.f, 1.f) * 255);
        data[int(4 * i + 3)] = (unsigned char) int((overlayAlpha.size() == overlay.size() ? overlayAlpha[i] : 1.f) * 255.f);       // Alpha
    }
    painter.drawImage(0, 0, QImage(data, image.sizeX, image.sizeY, QImage::Format_ARGB32));
    painter.end();
    return img;*/
    return this->computeDisplayedImage({{"", overlay}}, {{"", overlayAlpha}}, {{"", true}}, this->getImage().getDimensions());
}

QImage PlotImageData::computeDisplayedImage(const std::map<std::string, GridV3> &overlays, const std::map<std::string, GridF> &overlayAlphas, const std::map<std::string, bool>& displayedOverlays, const Vector3i &imgSize) const
{
    Vector3i largestDimensions = imgSize;
    for (auto& [name, over] : overlays) {
        largestDimensions.x() = std::max(largestDimensions.x(), (int)over.sizeX);
        largestDimensions.y() = std::max(largestDimensions.y(), (int)over.sizeY);
    }
    QImage img = this->computeDisplayedImage(largestDimensions);
    QPainter painter = QPainter(&img);

    for (auto& [name, over] : overlays) {
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

PlotVectorData *PlotVectorData::setField(const GridV3 &field)
{
    this->field = field;
    return this;
}

std::pair<GridV3, GridF> PlotVectorData::getFieldImageAndAlpha(const Vector3i &imgSize, const Vector3i& numberOfCells) const
{
    return PlotVectorData::createFieldImageAndAlpha(this->field, imgSize, numberOfCells);
}

std::pair<GridV3, GridF> PlotVectorData::createFieldImageAndAlpha(const GridV3 &field, Vector3i imgSize, const Vector3i &numberOfCells, const Vector3 &backgroundColor)
{
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
    Vector3 ratio = Vector3((float)imgSize.x() / (float)reducedSize.x(), (float)imgSize.y() / (float)reducedSize.y(), 1.f);
    reduced.iterateParallel([&] (const Vector3& _p) {
        Vector3 p = _p + Vector3(.5f, .5f);
        // AABBox cell((p - Vector3(.5f, .5f, 1)) * ratio, (p + Vector3(.5f, .5f, 1)) * ratio); // Added an depth (z) to avoid issue on the intersection computation
        Vector3 vec = reduced.interpolate(p);
        if (!vec.isValid()) return;
        float mag = vec.norm();
        if (mag < 1e-5) return;
        Vector3 dir = vec / mag;
        Vector3 color(1, 1, 1);
        float relativeMag = 1.f;
        if (std::abs(minMag - maxMag) > 1e-5) {
            relativeMag = interpolation::linear(mag, minMag, maxMag);
            color = colorPalette(relativeMag, {Vector3(70.f, 0.f, 100.f) / 255.f, Vector3(30.f, 160.f, 130.f) / 255.f, Vector3(255.f, 250.f, 0.f)/255.f});
        }
        bool valid = dir.xy().norm2() > 1e-5;
        if (!valid) return;

        Vector3 startLine = (p - dir * interpolation::inv_linear(relativeMag, .5f, 1.f)) * ratio;
        Vector3 endLine = (p + dir * interpolation::inv_linear(relativeMag, .5f, 1.f))  * ratio;
        float length = (endLine - startLine).norm();

        img = PlotVectorData::drawLine(img, color, startLine, endLine);
        alpha = PlotVectorData::drawLine(alpha, 1.f, startLine, endLine);

        img = PlotVectorData::drawLine(img, color, endLine, endLine - dir.rotated(deg2rad(20), Vector3(0, 0, 1)) * length * .3f);
        alpha = PlotVectorData::drawLine(alpha, 1.f, endLine, endLine - dir.rotated(deg2rad(20), Vector3(0, 0, 1)) * length * .3f);

        img = PlotVectorData::drawLine(img, color, endLine, endLine - dir.rotated(deg2rad(-20), Vector3(0, 0, 1)) * length * .3f);
        alpha = PlotVectorData::drawLine(alpha, 1.f, endLine, endLine - dir.rotated(deg2rad(-20), Vector3(0, 0, 1)) * length * .3f);
        /*float halfLength = (endLine - startLine).norm() * .5f;
        for (int i = 0; i < (int) (halfLength); i++) {
            const Vector3i coordA = p * ratio + dir * i;
            const Vector3i coordB = p * ratio - dir * i;
            img[coordA] = color;
            img[coordB] = color;
            alpha[coordA] = 1.f;
            alpha[coordB] = 1.f;
        }*/
        /*
        // Draw head
        const Vector3 dirHeadA = -dir.rotated(deg2rad(20), Vector3(0, 0, 1));
        const Vector3 dirHeadB = -dir.rotated(deg2rad(-20), Vector3(0, 0, 1));
        for (int i = 0; i < (int) (halfLength * 2.f * .33f); i++) {
            const Vector3i coord = endLine; //p * ratio + dir * maxI;
            img[coord + dirHeadA * i] = color;
            alpha[coord + dirHeadA * i] = 1.f;
            img[coord + dirHeadB * i] = color;
            alpha[coord + dirHeadB * i] = 1.f;
        }*/

    });
    return {img, alpha};
}
