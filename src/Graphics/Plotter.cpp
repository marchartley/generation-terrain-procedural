#include "Plotter.h"
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
    //    while (!this->chart()->axes().empty()) {
    //        this->chart()->removeAxis(this->chart()->axes().front());
    //    }

    //    for (auto& labels : this->graphicLabels)
    //        for (auto& lab : labels)
    //            delete lab;
    //    this->graphicLabels.clear();
    //    QTransform prevState = this->transform();
    //    auto prevState = this->chart()->transformations();
    //    this->chart()->removeAllSeries();
    //    while (!this->chart()->axes().empty()) {
    //        this->chart()->removeAxis(this->chart()->axes().front());
    //    }
    if (!title.empty())
        this->chart()->setTitle(QString::fromStdString(title));

    for (auto& labels : this->_dataModel->graphicLabels)
        for (auto& lab : labels)
            delete lab;
    this->_dataModel->graphicLabels.clear();

    if (!this->_dataModel->displayedImage.empty() && this->_dataModel->backImage && !this->_dataModel->backImage->isNull()) {
        int width = static_cast<int>(this->chart()->plotArea().width());
        int height = static_cast<int>(this->chart()->plotArea().height());
        int ViewW = static_cast<int>(this->width());
        int ViewH = static_cast<int>(this->height());

        //scale the image to fit plot area
        QImage scaledImage = this->_dataModel->backImage->scaled(QSize(width, height), Qt::IgnoreAspectRatio, Qt::TransformationMode::FastTransformation); // SmoothTransformation);
        //        *backImage = backImage->scaled(QSize(width, height));

        //We have to translate the image because setPlotAreaBackGround
        //starts the image in the top left corner of the view not the
        //plot area. So, to offset we will make a new image the size of
        //view and offset our image within that image with white
        QImage translated(ViewW, ViewH, QImage::Format_ARGB32);
        translated.fill(Qt::white);
        QPainter painter(&translated);
        QPointF TopLeft = this->chart()->plotArea().topLeft();
        painter.drawImage(TopLeft, scaledImage);

        //Display image in background
        //        this->chart()->setPlotAreaBackgroundBrush(scaledImage);
        this->chart()->setPlotAreaBackgroundBrush(translated);
        this->chart()->setPlotAreaBackgroundVisible(true);
    }
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
    //    this->setTransform(prevState);
    //    this->chart()->setTransformations(prevState);
    //    while (!this->chart()->axes().empty()) {
    //        this->chart()->removeAxis(this->chart()->axes().front());
    //    }

    this->chart()->createDefaultAxes();
    //    this->chart()->zoomOut();
    //    this->update();
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
        std::cout << "Removed a call to Plotter::draw() here..." << std::endl;
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
    } else {
        this->unlockView();
    }

    // if (!this->_dataModel->displayedImage.empty()) {
        // Q_EMIT clickedOnImage(pos * this->_dataModel->displayedImage.getDimensions(), this->_dataModel->displayedImage(pos * this->_dataModel->displayedImage.getDimensions()));
    // }
    return this;
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
    Q_EMIT this->clickedOnValue(vecLocal);

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
    this->selectData(Vector3::invalid());
    Q_EMIT this->clickedOnValue(Vector3(false)); // "unclick"
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

bool Chart::gestureEvent(QGestureEvent *event)
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

AbstractPlotter::AbstractPlotter(std::string name, QWidget *parent) : AbstractPlotter(name, new ChartView(new Chart()), parent)
{
}

AbstractPlotter::AbstractPlotter(std::string name, ChartView *chartView, QWidget *parent) : QDialog(parent), chartView(chartView), name(name)
{
    if (this->chartView == nullptr)
        this->chartView = new ChartView(new Chart());

    this->dataModel = new PlotModel();

    auto layout = new InterfaceUI(new QVBoxLayout());
    auto mainLayout = new InterfaceUI(new QHBoxLayout());

//    auto right = new QVBoxLayout();
    toolsInterface = new InterfaceUI(new QVBoxLayout(), "Tools");
    viewOptionsInterface = new InterfaceUI(new QVBoxLayout(), "View options");
    saveCopyInterface = new InterfaceUI(new QVBoxLayout(), "Save/Copy");
    infosInterface = new InterfaceUI(new QHBoxLayout(), "Infos");

    this->chartView->setRenderHint(QPainter::Antialiasing);
    this->chartView->chart()->legend()->setMarkerShape(QLegend::MarkerShapeFromSeries);
    auto viewAndCopyInterface = new InterfaceUI(new QVBoxLayout());
    viewAndCopyInterface->add(std::vector<UIElement*>({viewOptionsInterface, saveCopyInterface}));
//    this->chartView->setMaximumSize(10000, 10000);
//    this->chartView->chart()->setMaximumSize(10000, 10000);
    //    this->chartView->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    //    this->chartView->chart()->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    //    this->mouseInfoLabel->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Maximum);

    //    this->setWindowModality(Qt::WindowModality::NonModal);
    //    this->setModal(false);

    layout->add(std::vector<UIElement*>({mainLayout, infosInterface}));

    mainLayout->add({toolsInterface, new UIElement(this->chartView), viewAndCopyInterface});


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

    this->resize(800, 600);
    this->updateGeometry();

    this->setWindowTitle(QString::fromStdString(toCapitalize(this->name)));

    // QObject::connect(this->chartView, &ChartView::clickedOnValue, this, &AbstractPlotter::selectData);
    QObject::connect(this->chartView, &ChartView::mouseMoved, this, [&](const Vector3& pos, const Vector3& prevPos, QMouseEvent* e){
        this->displayInfoUnderMouse(pos);
        if (!this->dataModel->displayedImage.empty()) {
            Q_EMIT this->movedOnImage(pos * dataModel->displayedImage.getDimensions(), prevPos * dataModel->displayedImage.getDimensions(), e);
        }
    });
    QObject::connect(this->chartView->chart(), &QChart::geometryChanged, this, &AbstractPlotter::draw);
    QObject::connect(this->chartView, &ChartView::clickedOnValue, this, [&](const Vector3& pos) {
        if (this->dataModel->displayedImage.size() > 0)
            Q_EMIT this->clickedOnImage(pos, this->dataModel->displayedImage.at(pos * this->dataModel->displayedImage.getDimensions()));
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

AbstractPlotter* AbstractPlotter::addImage(GridV3 image)
{
    this->dataModel->addImage(image);
    return this;
}

AbstractPlotter* AbstractPlotter::addImage(const GridF &image)
{
    GridV3 copy(image.getDimensions());
    for (size_t i = 0; i < copy.size(); i++) {
        float val = image[i];
        copy[i] = Vector3(val, val, val);
    }
    return this->addImage(copy);
}

AbstractPlotter* AbstractPlotter::addImage(const Matrix3<double> &image)
{
    return this->addImage((GridF)image);
}

AbstractPlotter* AbstractPlotter::addImage(const GridI &image)
{
    return this->addImage((GridF)image);
}

GridV3 AbstractPlotter::computeVectorFieldRendering(const GridV3 &field, float reductionFactor, Vector3 imgSize) const
{
    if (!imgSize.isValid())
        imgSize = field.getDimensions();
    imgSize.z() = 1;
    GridV3 img(imgSize);
    Vector3 reducedSize = (field.getDimensions() * reductionFactor).roundedUp(); //(30, 30, 1);
    reducedSize.z() = 1;
    Vector3 ratio = imgSize / reducedSize;
    GridV3 reduced = field.resize(reducedSize);

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

    reduced.iterateParallel([&] (const Vector3& p) {
        AABBox cell(p - Vector3(.45, .45, 0) * ratio, p + Vector3(.45, .45, 0) * ratio);
        float mag = reduced(p).norm();
        if (mag < 1e-5) return;
        Vector3 dir = reduced(p) / mag;
        Vector3 color(1, 1, 1);
        if (std::abs(minMag - maxMag) > 1e-5) {
            float relativeMag = interpolation::linear(mag, minMag, maxMag);
            color = colorPalette(relativeMag, {Vector3(0, 0, 1), Vector3(1, 1, 1), Vector3(1, 0, 0)});
        }
        bool valid = true;
        for (int i = 0; valid; i++) {
            valid = false;
            if (cell.containsXY(p + dir * i)) {
                img((p + Vector3(.5, .5)) * ratio + dir * i) = color;
                valid = true;
            }
        }
    });
    return img;
}

AbstractPlotter* AbstractPlotter::addVectorField(const GridV3 &field, float reductionFactor, Vector3 imgSize, float opacity)
{
    GridV3 img = this->computeVectorFieldRendering(field, reductionFactor, imgSize);
    if (this->dataModel->displayedImage.size() > 0) {
        img = this->dataModel->displayedImage.resize(img.getDimensions()) + img * (opacity);
    }
    return this->addImage(img);
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
    if (this->dataModel->displayedImage.size() > 0) {
        img = this->dataModel->displayedImage.resize(img.getDimensions()) + img * (opacity);
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
    this->draw();
    QDialog::show();
    return this;
}

AbstractPlotter* AbstractPlotter::updateUI()
{
    blockSignals(true);

    this->toolsInterface->update();
    this->viewOptionsInterface->update();
    this->saveCopyInterface->update();
    this->infosInterface->update();
    blockSignals(false);
    return this;
}

int AbstractPlotter::exec()
{
    this->draw();

    return QDialog::exec();

}

AbstractPlotter* AbstractPlotter::saveFig(std::string filename)
{
    QPixmap p = this->chartView->grab();
    if (this->dataModel->backImage)
        p = QPixmap::fromImage(*this->dataModel->backImage);
    p.save(QString::fromStdString(filename), "PNG");
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
    // this->chartView->reset();
    /*
    this->chartView->chart()->removeAllSeries();
    while (!this->chartView->chart()->axes().empty()) {
        this->chartView->chart()->removeAxis(this->chartView->chart()->axes().front());
    }
    if (!title.empty())
        this->chartView->chart()->setTitle(QString::fromStdString(title));

    for (auto& labels : this->graphicLabels)
        for (auto& lab : labels)
            delete lab;
    this->graphicLabels.clear();

//    QPushButton* saveButton;
//    ChartView* chartView;
    backImage = nullptr;
    title = "";
    plot_data.clear();
    plot_names.clear();
    plot_colors.clear();
    scatter_data.clear();
    scatter_labels.clear();
    scatter_colors.clear();
    scatter_names.clear();
    graphicLabels.clear();

    selectedScatterData.clear();
    selectedPlotData.clear();
    */
    return this;
}
/*
AbstractPlotter* AbstractPlotter::updateLabelsPositions()
{
//    this->blockSignals(true);
    if (!this->selectedPlotData.empty() || !this->selectedScatterData.empty()) {
        QPointF qNewPoint = this->chartView->chart()->mapToValue(this->chartView->previousMousePos);
        Vector3 newPoint = Vector3(qNewPoint.x()(), qNewPoint.y()());
        for (auto& [iPlot, iPoint] : this->selectedPlotData)
            this->plot_data[iPlot][iPoint] = newPoint;
        for (auto& [iPlot, iPoint] : this->selectedScatterData)
            this->scatter_data[iPlot][iPoint] = newPoint;
    }

    for (size_t iScatter = 0; iScatter < this->scatter_labels.size(); iScatter++) {
        for (size_t iPoint = 0; iPoint < this->scatter_labels[iScatter].size(); iPoint++) {
//                this->graphicLabels[iScatter][iPoint]->setPos(QPointF(this->scatter_data[iScatter][iPoint].first, this->scatter_data[iScatter][iPoint].second)); // this->chartView->chart()->mapToPosition(QPointF(this->scatter_data[iScatter][iPoint].first, this->scatter_data[iScatter][iPoint].second)));
            this->graphicLabels[iScatter][iPoint]->setPos(this->chartView->chart()->mapToPosition(QPointF(this->scatter_data[iScatter][iPoint].x(), this->scatter_data[iScatter][iPoint].y())));
        }
    }
    if (!this->selectedPlotData.empty() || !this->selectedScatterData.empty()) {
        this->draw();
    }
//    this->blockSignals(false);
    return this;
}

AbstractPlotter* AbstractPlotter::selectData(const Vector3& pos)
{
    if (!pos.isValid()) return this;

    float minDist = 0.05f;
    this->selectedScatterData.clear();
    this->selectedPlotData.clear();

    if (pos.isValid()) {
        for (size_t i = 0; i < this->plot_data.size(); i++) {
            for (size_t j = 0; j < this->plot_data[i].size(); j++) {
                if ((plot_data[i][j] - pos).norm2() < minDist*minDist)
                    this->selectedPlotData.push_back({i, j});
            }
        }
        for (size_t i = 0; i < this->scatter_data.size(); i++) {
            for (size_t j = 0; j < this->scatter_data[i].size(); j++) {
                if ((scatter_data[i][j] - pos).norm2() < minDist*minDist)
                    this->selectedScatterData.push_back({i, j});
            }
        }
    }

    if (!this->selectedPlotData.empty() || !this->selectedScatterData.empty()) {
        this->chartView->lockView();
    } else {
        this->chartView->unlockView();
    }

    if (!this->displayedImage.empty()) {
        Q_EMIT clickedOnImage(pos * displayedImage.getDimensions(), this->displayedImage(pos * displayedImage.getDimensions()));
    }
    return this;
}
*/

AbstractPlotter* AbstractPlotter::displayInfoUnderMouse(const Vector3 &relativeMousePos)
{
    if (this->dataModel->displayedImage.empty() || relativeMousePos.minComp() < 0.f || relativeMousePos.maxComp() > 1.f)
        return this;
    std::ostringstream oss;
    Vector3 size = this->dataModel->displayedImage.getDimensions();
    Vector3 position = relativeMousePos * size;
    Vector3 value = this->dataModel->displayedImage(position);
    oss << "Mouse pos: " << int(position.x()) << ", " << int(position.y()) << " -- Value : (" << value.x() << ", " << value.y() << ", " << value.z() << ") ";
    this->mouseInfoLabel->setText(QString::fromStdString(oss.str()));
    return this;
}













Plotter::Plotter(std::string name, QWidget *parent) : Plotter(name, new ChartView(new Chart()), parent)
{
}

Plotter::Plotter(std::string name, ChartView *chartView, QWidget *parent) : AbstractPlotter(name, chartView, parent)
{
    /*if (this->chartView == nullptr)
        this->chartView = new ChartView(new Chart());

    auto layout = new QHBoxLayout();
    auto left = new QVBoxLayout();
    //    auto right = new QVBoxLayout();
    interfaceButtons = new InterfaceUI(new QVBoxLayout(), "Display");

    this->chartView->setRenderHint(QPainter::Antialiasing);
    this->chartView->chart()->legend()->setMarkerShape(QLegend::MarkerShapeFromSeries);
    //    this->chartView->setMaximumSize(10000, 10000);
    //    this->chartView->chart()->setMaximumSize(10000, 10000);

    //    this->saveButton = new ButtonElement("Save");
    auto saveButton = new ButtonElement("Save");
    auto copyToClipboardButton = new ButtonElement("Copy");
    this->mouseInfoLabel = new QLabel("");


    //    this->chartView->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    //    this->chartView->chart()->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    //    this->mouseInfoLabel->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Maximum);

    */

    auto normalizeModeButton = new CheckboxElement("Normalize", normalizedMode);
    auto absoluteModeButton = new CheckboxElement("Absolute", absoluteMode);
    this->rangeValuesWidget = new RangeSliderElement("Values", -1000, 1000, 0.01f, minValueToDisplay, maxValueToDisplay, Qt::Vertical);
    auto rangeActiveCheckbox = new CheckboxElement("Filter", clampValues);

    absoluteModeButton->setOnChecked([&](bool toggled) {
        this->addImage(this->dataModel->displayedImage);
        this->draw();
    });
    normalizeModeButton->setOnChecked([&](bool toggled) {
        this->addImage(this->dataModel->displayedImage);
        this->draw();
    });
    rangeValuesWidget->setOnValueChanged([&](float minVal, float maxVal) {
        this->addImage(this->dataModel->displayedImage);
        this->draw();
    });

    CheckboxElement* displayRButton = new CheckboxElement("R", this->displayR);
    CheckboxElement* displayGButton = new CheckboxElement("G", this->displayG);
    CheckboxElement* displayBButton = new CheckboxElement("B", this->displayB);

    displayRButton->setOnChecked([&](bool toggled) { this->addImage(this->dataModel->displayedImage); this->draw(); });
    displayGButton->setOnChecked([&](bool toggled) { this->addImage(this->dataModel->displayedImage); this->draw(); });
    displayBButton->setOnChecked([&](bool toggled) { this->addImage(this->dataModel->displayedImage); this->draw(); });

    this->viewOptionsInterface->add(std::vector<UIElement*>{normalizeModeButton, absoluteModeButton, rangeValuesWidget, displayRButton, displayGButton, displayBButton});
}

Plotter *Plotter::getInstance(std::string name)
{
    if (name == "") name = Plotter::defaultName;
    if (Plotter::instances.count(name) == 0) {
        //        std::cerr << "Plotter has not been initialized with function Plotter::init()" << std::endl;
        Plotter::instances[name] = Plotter::init(name);
    }
    return dynamic_cast<Plotter*>(Plotter::instances[name]);
}

Plotter *Plotter::init(std::string name, ChartView *chartView, QWidget *parent)
{
    if (Plotter::instances.count(name))
        delete Plotter::instances[name];
    Plotter::instances[name] = new Plotter(name, chartView, parent);
    return Plotter::getInstance(name);
}
/*
Plotter* Plotter::draw()
{
    this->chartView->chart()->removeAllSeries();
    //    while (!this->chartView->chart()->axes().empty()) {
    //        this->chartView->chart()->removeAxis(this->chartView->chart()->axes().front());
    //    }

    //    for (auto& labels : this->graphicLabels)
    //        for (auto& lab : labels)
    //            delete lab;
    //    this->graphicLabels.clear();
    //    QTransform prevState = this->chartView->transform();
    //    auto prevState = this->chartView->chart()->transformations();
    //    this->chartView->chart()->removeAllSeries();
    //    while (!this->chartView->chart()->axes().empty()) {
    //        this->chartView->chart()->removeAxis(this->chartView->chart()->axes().front());
    //    }
    if (!title.empty())
        this->chartView->chart()->setTitle(QString::fromStdString(title));

    for (auto& labels : this->graphicLabels)
        for (auto& lab : labels)
            delete lab;
    this->graphicLabels.clear();

    if (!this->displayedImage.empty() && this->backImage && !this->backImage->isNull()) {
        int width = static_cast<int>(this->chartView->chart()->plotArea().width());
        int height = static_cast<int>(this->chartView->chart()->plotArea().height());
        int ViewW = static_cast<int>(chartView->width());
        int ViewH = static_cast<int>(chartView->height());

        //scale the image to fit plot area
        QImage scaledImage = backImage->scaled(QSize(width, height), Qt::IgnoreAspectRatio, Qt::TransformationMode::FastTransformation); // SmoothTransformation);
        //        *backImage = backImage->scaled(QSize(width, height));

        //We have to translate the image because setPlotAreaBackGround
        //starts the image in the top left corner of the view not the
        //plot area. So, to offset we will make a new image the size of
        //view and offset our image within that image with white
        QImage translated(ViewW, ViewH, QImage::Format_ARGB32);
        translated.fill(Qt::white);
        QPainter painter(&translated);
        QPointF TopLeft = this->chartView->chart()->plotArea().topLeft();
        painter.drawImage(TopLeft, scaledImage);

        //Display image in background
        //        this->chartView->chart()->setPlotAreaBackgroundBrush(scaledImage);
        this->chartView->chart()->setPlotAreaBackgroundBrush(translated);
        this->chartView->chart()->setPlotAreaBackgroundVisible(true);
    }
    for (size_t i = 0; i < this->plot_data.size(); i++) {
        QLineSeries *series = new QLineSeries();
        if (this->plot_names.size() > 0 && this->plot_names.size() == this->plot_data.size())
            series->setName(QString::fromStdString(this->plot_names[i]));
        for (auto pos : this->plot_data[i])
            series->append(pos.x, pos.y);
        //        if (PlotColorToQColor.count(this->plot_colors[i]))
        //            series->setColor(PlotColorToQColor.at(this->plot_colors[i]));
        series->setColor(this->plot_colors[i]);
        this->chartView->chart()->addSeries(series);
        if (series->name().isEmpty()) {
            this->chartView->chart()->legend()->markers(series)[0]->setVisible(false);
        }
    }
    for (size_t i = 0; i < this->scatter_data.size(); i++) {
        QScatterSeries *series = new QScatterSeries();
        if (this->scatter_names.size() > 0 && this->scatter_names.size() == this->scatter_data.size())
            series->setName(QString::fromStdString(this->scatter_names[i]));
        for (size_t j = 0; j < this->scatter_data[i].size(); j++) { //(auto pos : this->scatter_data[i]) {
            auto pos = this->scatter_data[i][j];
            series->append(pos.x, pos.y);
        }
        this->chartView->chart()->addSeries(series);

        if (series->name().isEmpty()) {
            this->chartView->chart()->legend()->markers(series)[0]->setVisible(false);
        }
    }
    for (size_t iScatter = 0; iScatter < this->scatter_labels.size(); iScatter++) {
        this->graphicLabels.push_back(std::vector<QGraphicsTextItem*>());
        for (size_t iPoint = 0; iPoint < this->scatter_labels[iScatter].size(); iPoint++) {
            QGraphicsTextItem *itm = new QGraphicsTextItem(QString::fromStdString(this->scatter_labels[iScatter][iPoint]), this->chartView->chart());
            this->graphicLabels[iScatter].push_back(itm);
        }
    }
    //    this->chartView->setTransform(prevState);
    //    this->chartView->chart()->setTransformations(prevState);
    //    while (!this->chartView->chart()->axes().empty()) {
    //        this->chartView->chart()->removeAxis(this->chartView->chart()->axes().front());
    //    }

    this->chartView->chart()->createDefaultAxes();
    //    this->chartView->chart()->zoomOut();
    //    this->chartView->update();
    return this;
}
*/
Plotter* Plotter::updateUI()
{
    AbstractPlotter::updateUI();
    return this;
}

Plotter* Plotter::setNormalizedModeImage(bool normalize)
{
    this->normalizedMode = normalize;
    this->addImage(this->dataModel->displayedImage);
    return updateUI();
}

Plotter* Plotter::setAbsoluteModeImage(bool absolute)
{
    this->absoluteMode = absolute;
    this->addImage(this->dataModel->displayedImage);
    return updateUI();
}

Plotter *Plotter::setFilteredValuesImage(bool filtered)
{
    this->clampValues = filtered;
    this->addImage(this->dataModel->displayedImage);
    return updateUI();
}
/*
Plotter* Plotter::reset()
{
    AbstractPlotter::reset();

    return this;
}
*/








/*
std::map<std::string, Plotter*> Plotter::instances = std::map<std::string, Plotter*>();
std::string Plotter::defaultName = "default";

Plotter::Plotter(std::string name, QWidget *parent) : Plotter(name, new ChartView(new Chart()), parent)
{
}

Plotter::Plotter(std::string name, ChartView *chartView, QWidget *parent) : QDialog(parent), chartView(chartView), name(name)
{
    if (this->chartView == nullptr)
        this->chartView = new ChartView(new Chart());

    auto layout = new QHBoxLayout();
    auto left = new QVBoxLayout();
    //    auto right = new QVBoxLayout();
    interfaceButtons = new InterfaceUI(new QVBoxLayout(), "Display");

    this->chartView->setRenderHint(QPainter::Antialiasing);
    this->chartView->chart()->legend()->setMarkerShape(QLegend::MarkerShapeFromSeries);
    //    this->chartView->setMaximumSize(10000, 10000);
    //    this->chartView->chart()->setMaximumSize(10000, 10000);

    auto normalizeModeButton = new CheckboxElement("Normalize", normalizedMode);
    auto absoluteModeButton = new CheckboxElement("Absolute", absoluteMode);
    this->rangeValuesWidget = new RangeSliderElement("Values", -1000, 1000, 0.01f, minValueToDisplay, maxValueToDisplay, Qt::Vertical);
    auto rangeActiveCheckbox = new CheckboxElement("Filter", clampValues);
    //    this->saveButton = new ButtonElement("Save");
    auto saveButton = new ButtonElement("Save");
    auto copyToClipboardButton = new ButtonElement("Copy");
    this->mouseInfoLabel = new QLabel("");


    //    this->chartView->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    //    this->chartView->chart()->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    //    this->mouseInfoLabel->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Maximum);


    absoluteModeButton->setOnChecked([&](bool toggled) {
        this->addImage(displayedImage);
        this->draw();
    });
    normalizeModeButton->setOnChecked([&](bool toggled) {
        this->addImage(displayedImage);
        this->draw();
    });
    saveButton->setOnClick([&]() {
        QString q_filename = QFileDialog::getSaveFileName(this, QString("Save plot"));
        if (!q_filename.isEmpty())
            saveFig(q_filename.toStdString());
    });
    copyToClipboardButton->setOnClick([&]() { this->copyToClipboard(); });

    rangeValuesWidget->setOnValueChanged([&](float minVal, float maxVal) {
        this->addImage(displayedImage);
        this->draw();
    });

    CheckboxElement* displayRButton = new CheckboxElement("R", this->displayR);
    CheckboxElement* displayGButton = new CheckboxElement("G", this->displayG);
    CheckboxElement* displayBButton = new CheckboxElement("B", this->displayB);

    displayRButton->setOnChecked([&](bool toggled) { this->addImage(displayedImage); this->draw(); });
    displayGButton->setOnChecked([&](bool toggled) { this->addImage(displayedImage); this->draw(); });
    displayBButton->setOnChecked([&](bool toggled) { this->addImage(displayedImage); this->draw(); });

    InterfaceUI* RGBSelection = new InterfaceUI(new QHBoxLayout);
    RGBSelection->add({displayRButton, displayGButton, displayBButton});
    RGBSelection->get()->layout()->setSpacing(0);
    RGBSelection->get()->layout()->setContentsMargins(0, 0, 0, 0);

    left->addWidget(RGBSelection->get());
    left->addWidget(this->chartView);
    left->addWidget(this->mouseInfoLabel);
    interfaceButtons->add(normalizeModeButton);
    interfaceButtons->add(absoluteModeButton);
    interfaceButtons->add(rangeActiveCheckbox);
    interfaceButtons->add(this->rangeValuesWidget);
    interfaceButtons->add(saveButton);
    interfaceButtons->add(copyToClipboardButton);

    //    this->setWindowModality(Qt::WindowModality::NonModal);
    //    this->setModal(false);

    layout->addItem(left);
    layout->addWidget(interfaceButtons->get());
    this->setLayout(layout);

    this->resize(800, 600);
    this->updateGeometry();

    this->setWindowTitle(QString::fromStdString(toCapitalize(this->name)));

    QObject::connect(this->chartView, &ChartView::clickedOnValue, this, &Plotter::selectData);
    QObject::connect(this->chartView->chart(), &QChart::geometryChanged, this, &Plotter::updateLabelsPositions);
    QObject::connect(this->chartView->chart(), &QChart::plotAreaChanged, this, &Plotter::updateLabelsPositions);
    QObject::connect(this->chartView, &ChartView::updated, this, &Plotter::updateLabelsPositions);
    QObject::connect(this->chartView, &ChartView::mouseMoved, this, [&](const Vector3& pos, const Vector3& prevPos, QMouseEvent* e){
        this->displayInfoUnderMouse(pos);
        if (!this->displayedImage.empty()) {
            Q_EMIT this->movedOnImage(pos * displayedImage.getDimensions(), prevPos * displayedImage.getDimensions(), e);
        }
    });
    QObject::connect(this->chartView->chart(), &QChart::geometryChanged, this, &Plotter::draw);
}

Plotter *Plotter::getInstance(std::string name)
{
    if (name == "") name = Plotter::defaultName;
    if (Plotter::instances.count(name) == 0) {
        //        std::cerr << "Plotter has not been initialized with function Plotter::init()" << std::endl;
        Plotter::instances[name] = Plotter::init(name);
    }
    return Plotter::instances[name];
}

Plotter *Plotter::init(std::string name, ChartView *chartView, QWidget *parent)
{
    if (Plotter::instances.count(name))
        delete Plotter::instances[name];
    Plotter::instances[name] = new Plotter(name, chartView, parent);
    return Plotter::getInstance(name);
}

Plotter* Plotter::addPlot(std::vector<float> data, std::string name, QColor color)
{
    std::vector<Vector3> _data;
    for (unsigned int i = 0; i < data.size(); i++) {
        _data.push_back(Vector3(i, data[i]));
    }
    this->addPlot(_data, name, color);
    return this;
}

Plotter* Plotter::addPlot(std::vector<Vector3> data, std::string name, QColor color)
{
    this->plot_data.push_back(data);
    this->plot_names.push_back(name);
    this->plot_colors.push_back(color);
    return this;
}

Plotter *Plotter::addPlot(const BSpline &data, std::string name, QColor color)
{
    return this->addPlot(data.points, name, color);
}

Plotter* Plotter::addScatter(std::vector<float> data, std::string name, std::vector<std::string> labels, std::vector<QColor> colors)
{
    std::vector<Vector3> _data;
    for (unsigned int i = 0; i < data.size(); i++) {
        _data.push_back(Vector3(i, data[i]));
    }
    this->addScatter(_data, name, labels, colors);
    return this;
}

Plotter* Plotter::addScatter(std::vector<Vector3> data, std::string name, std::vector<std::string> labels, std::vector<QColor> colors)
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

Plotter* Plotter::addImage(GridV3 image)
{
    this->displayedImage = image;
    //    image = image.flip(false, true);
    if (image.empty()) return this;
    if (this->clampValues) {
        float min = std::numeric_limits<float>::max();
        float max = std::numeric_limits<float>::lowest();
        image.iterate([&](size_t i) {
            min = std::min(min, image[i].minComp());
            max = std::max(max, image[i].maxComp());
        });
        this->rangeValuesWidget->slider()->setMinimalValue(min);
        this->rangeValuesWidget->slider()->setMaximalValue(max);
        if (minValueToDisplay < min) minValueToDisplay = min;
        if (maxValueToDisplay > max) maxValueToDisplay = max;
        image.iterateParallel([&](size_t i) {
            for (int c = 0; c < 3; c++) {
                image[i][c] = std::clamp(image[i][c], this->minValueToDisplay, this->maxValueToDisplay);
            }
        });
    }
    if (this->absoluteMode) {
        image = image.abs();
    }
    if (this->normalizedMode) {
        for (int c = 0; c < 3; c++) {
            float min = std::numeric_limits<float>::max();
            float max = std::numeric_limits<float>::lowest();
            image.iterate([&](size_t i) {
                min = std::min(image[i][c], min);
                max = std::max(image[i][c], max);
            });
            float d = max - min;
            if (d == 0) {
                image.iterateParallel([&](size_t i) {
                    image[i][c] = 0.f;
                });
            } else {
                image.iterateParallel([&](size_t i) {
                    image[i][c] = (image[i][c] - min) / d;
                });
            }
        }
        //        image.normalize();
    }
    Vector3 colorFilter = Vector3((displayR ? 1.f : 0.f), (displayG ? 1.f : 0.f), (displayB ? 1.f : 0.f));
    image.iterateParallel([&](size_t i) {
        image[i] *= colorFilter;
    });
    unsigned char* data = new unsigned char[image.size() * 4];

    for (size_t i = 0; i < image.size(); ++i) {
        data[int(4 * i + 2)] = (unsigned char)(std::clamp(image[i].x, 0.f, 1.f) * 255);
        data[int(4 * i + 1)] = (unsigned char)(std::clamp(image[i].y, 0.f, 1.f) * 255);
        data[int(4 * i + 0)] = (unsigned char)(std::clamp(image[i].z, 0.f, 1.f) * 255);
        data[int(4 * i + 3)] = (unsigned char) 255;       // Alpha
    }

    if (this->backImage) {
        delete this->backImage;
    }
    this->backImage = new QImage(data, image.sizeX, image.sizeY, QImage::Format_ARGB32);
    //    *(this->backImage) = this->backImage->mirrored();
    return this;
}

Plotter* Plotter::addImage(const GridF &image)
{
    GridV3 copy(image.getDimensions());
    for (size_t i = 0; i < copy.size(); i++) {
        float val = image[i];
        copy[i] = Vector3(val, val, val);
    }
    return this->addImage(copy);
}

Plotter* Plotter::addImage(const Matrix3<double> &image)
{
    return this->addImage((GridF)image);
}

Plotter* Plotter::addImage(const GridI &image)
{
    return this->addImage((GridF)image);
}

GridV3 Plotter::computeVectorFieldRendering(const GridV3 &field, float reductionFactor, Vector3 imgSize) const
{
    if (!imgSize.isValid())
        imgSize = field.getDimensions();
    imgSize.z() = 1;
    GridV3 img(imgSize);
    Vector3 reducedSize = (field.getDimensions() * reductionFactor).roundedUp(); //(30, 30, 1);
    reducedSize.z() = 1;
    Vector3 ratio = imgSize / reducedSize;
    GridV3 reduced = field.resize(reducedSize);

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

    reduced.iterateParallel([&] (const Vector3& p) {
        AABBox cell(p - Vector3(.45, .45, 0) * ratio, p + Vector3(.45, .45, 0) * ratio);
        float mag = reduced(p).norm();
        if (mag < 1e-5) return;
        Vector3 dir = reduced(p) / mag;
        Vector3 color(1, 1, 1);
        if (std::abs(minMag - maxMag) > 1e-5) {
            float relativeMag = interpolation::linear(mag, minMag, maxMag);
            color = colorPalette(relativeMag, {Vector3(0, 0, 1), Vector3(1, 1, 1), Vector3(1, 0, 0)});
        }
        bool valid = true;
        for (int i = 0; valid; i++) {
            valid = false;
            if (cell.containsXY(p + dir * i)) {
                img((p + Vector3(.5, .5)) * ratio + dir * i) = color;
                valid = true;
            }
        }
    });
    return img;
}

Plotter* Plotter::addVectorField(const GridV3 &field, float reductionFactor, Vector3 imgSize, float opacity)
{
    GridV3 img = this->computeVectorFieldRendering(field, reductionFactor, imgSize);
    if (this->displayedImage.size() > 0) {
        img = displayedImage.resize(img.getDimensions()) + img * (opacity);
    }
    return this->addImage(img);
}

GridV3 Plotter::computeStreamLinesRendering(const GridV3 &field, Vector3 imgSize) const
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

Plotter* Plotter::addStreamLines(const GridV3 &field, Vector3 imgSize, float opacity)
{
    GridV3 img = computeStreamLinesRendering(field, imgSize);
    if (this->displayedImage.size() > 0) {
        img = displayedImage.resize(img.getDimensions()) + img * (opacity);
    }
    return this->addImage(img);
}

Plotter* Plotter::draw()
{
    this->chartView->chart()->removeAllSeries();
    //    while (!this->chartView->chart()->axes().empty()) {
    //        this->chartView->chart()->removeAxis(this->chartView->chart()->axes().front());
    //    }

    //    for (auto& labels : this->graphicLabels)
    //        for (auto& lab : labels)
    //            delete lab;
    //    this->graphicLabels.clear();
    //    QTransform prevState = this->chartView->transform();
    //    auto prevState = this->chartView->chart()->transformations();
    //    this->chartView->chart()->removeAllSeries();
    //    while (!this->chartView->chart()->axes().empty()) {
    //        this->chartView->chart()->removeAxis(this->chartView->chart()->axes().front());
    //    }
    if (!title.empty())
        this->chartView->chart()->setTitle(QString::fromStdString(title));

    for (auto& labels : this->graphicLabels)
        for (auto& lab : labels)
            delete lab;
    this->graphicLabels.clear();

    if (!this->displayedImage.empty() && this->backImage && !this->backImage->isNull()) {
        int width = static_cast<int>(this->chartView->chart()->plotArea().width());
        int height = static_cast<int>(this->chartView->chart()->plotArea().height());
        int ViewW = static_cast<int>(chartView->width());
        int ViewH = static_cast<int>(chartView->height());

        //scale the image to fit plot area
        QImage scaledImage = backImage->scaled(QSize(width, height), Qt::IgnoreAspectRatio, Qt::TransformationMode::FastTransformation); // SmoothTransformation);
        //        *backImage = backImage->scaled(QSize(width, height));

        //We have to translate the image because setPlotAreaBackGround
        //starts the image in the top left corner of the view not the
        //plot area. So, to offset we will make a new image the size of
        //view and offset our image within that image with white
        QImage translated(ViewW, ViewH, QImage::Format_ARGB32);
        translated.fill(Qt::white);
        QPainter painter(&translated);
    Plotter* updateLabelsPositions();
    Plotter* selectData(const Vector3& pos);
    Plotter* displayInfoUnderMouse(const Vector3& relativeMousePos);
    Plotter* draw();
    Plotter* show();
        QPointF TopLeft = this->chartView->chart()->plotArea().topLeft();
        painter.drawImage(TopLeft, scaledImage);

        //Display image in background
        //        this->chartView->chart()->setPlotAreaBackgroundBrush(scaledImage);
        this->chartView->chart()->setPlotAreaBackgroundBrush(translated);
        this->chartView->chart()->setPlotAreaBackgroundVisible(true);
    }
    for (size_t i = 0; i < this->plot_data.size(); i++) {
        QLineSeries *series = new QLineSeries();
        if (this->plot_names.size() > 0 && this->plot_names.size() == this->plot_data.size())
            series->setName(QString::fromStdString(this->plot_names[i]));
        for (auto pos : this->plot_data[i])
            series->append(pos.x, pos.y);
        //        if (PlotColorToQColor.count(this->plot_colors[i]))
        //            series->setColor(PlotColorToQColor.at(this->plot_colors[i]));
        series->setColor(this->plot_colors[i]);
        this->chartView->chart()->addSeries(series);
        if (series->name().isEmpty()) {
            this->chartView->chart()->legend()->markers(series)[0]->setVisible(false);
        }
    }
    for (size_t i = 0; i < this->scatter_data.size(); i++) {
        QScatterSeries *series = new QScatterSeries();
        if (this->scatter_names.size() > 0 && this->scatter_names.size() == this->scatter_data.size())
            series->setName(QString::fromStdString(this->scatter_names[i]));
        for (size_t j = 0; j < this->scatter_data[i].size(); j++) { //(auto pos : this->scatter_data[i]) {
            auto pos = this->scatter_data[i][j];
            series->append(pos.x, pos.y);
        }
        this->chartView->chart()->addSeries(series);

        if (series->name().isEmpty()) {
            this->chartView->chart()->legend()->markers(series)[0]->setVisible(false);
        }
    }
    for (size_t iScatter = 0; iScatter < this->scatter_labels.size(); iScatter++) {
        this->graphicLabels.push_back(std::vector<QGraphicsTextItem*>());
        for (size_t iPoint = 0; iPoint < this->scatter_labels[iScatter].size(); iPoint++) {
            QGraphicsTextItem *itm = new QGraphicsTextItem(QString::fromStdString(this->scatter_labels[iScatter][iPoint]), this->chartView->chart());
            this->graphicLabels[iScatter].push_back(itm);
        }
    }
    //    this->chartView->setTransform(prevState);
    //    this->chartView->chart()->setTransformations(prevState);
    //    while (!this->chartView->chart()->axes().empty()) {
    //        this->chartView->chart()->removeAxis(this->chartView->chart()->axes().front());
    //    }

    this->chartView->chart()->createDefaultAxes();
    //    this->chartView->chart()->zoomOut();
    //    this->chartView->update();
    return this;
}

Plotter* Plotter::show()
{
    this->draw();
    QDialog::show();
    return this;
}

Plotter* Plotter::updateUI()
{
    blockSignals(true);
    this->interfaceButtons->update();
    blockSignals(false);
    return this;
}

Plotter* Plotter::setNormalizedModeImage(bool normalize)
{
    this->normalizedMode = normalize;
    this->addImage(this->displayedImage);
    return updateUI();
}

Plotter* Plotter::setAbsoluteModeImage(bool absolute)
{
    this->absoluteMode = absolute;
    this->addImage(this->displayedImage);
    return updateUI();
}

Plotter *Plotter::setFilteredValuesImage(bool filtered)
{
    this->clampValues = filtered;
    this->addImage(this->displayedImage);
    return updateUI();
}

int Plotter::exec()
{
    this->draw();

    return QDialog::exec();

}

Plotter* Plotter::saveFig(std::string filename)
{
    QPixmap p = this->chartView->grab();
    if (this->backImage)
        p = QPixmap::fromImage(*backImage);
    p.save(QString::fromStdString(filename), "PNG");
    return this;
}

Plotter* Plotter::copyToClipboard()
{
    QPixmap p = this->chartView->grab();
    QApplication::clipboard()->setPixmap(p, QClipboard::Clipboard);
    return this;
}

void Plotter::resizeEvent(QResizeEvent *event)
{
    QDialog::resizeEvent(event);
    this->draw();
}

void Plotter::showEvent(QShowEvent *event)
{
    QDialog::showEvent(event);
    this->draw();
}

QTimer *Plotter::animate(std::function<void ()> callback, int interval_ms)
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

Plotter* Plotter::reset()
{
    auto prevState = this->chartView->chart()->transformations();
    this->chartView->chart()->removeAllSeries();
    while (!this->chartView->chart()->axes().empty()) {
        this->chartView->chart()->removeAxis(this->chartView->chart()->axes().front());
    }
    if (!title.empty())
        this->chartView->chart()->setTitle(QString::fromStdString(title));

    for (auto& labels : this->graphicLabels)
        for (auto& lab : labels)
            delete lab;
    this->graphicLabels.clear();

    //    QPushButton* saveButton;
    //    ChartView* chartView;
    backImage = nullptr;
    title = "";
    plot_data.clear();
    plot_names.clear();
    plot_colors.clear();
    scatter_data.clear();
    scatter_labels.clear();
    scatter_colors.clear();
    scatter_names.clear();
    graphicLabels.clear();

    selectedScatterData.clear();
    selectedPlotData.clear();

    return this;
}

Plotter* Plotter::updateLabelsPositions()
{
    //    this->blockSignals(true);
    if (!this->selectedPlotData.empty() || !this->selectedScatterData.empty()) {
        QPointF qNewPoint = this->chartView->chart()->mapToValue(this->chartView->previousMousePos);
        Vector3 newPoint = Vector3(qNewPoint.x(), qNewPoint.y());
        for (auto& [iPlot, iPoint] : this->selectedPlotData)
            this->plot_data[iPlot][iPoint] = newPoint;
        for (auto& [iPlot, iPoint] : this->selectedScatterData)
            this->scatter_data[iPlot][iPoint] = newPoint;
    }

    for (size_t iScatter = 0; iScatter < this->scatter_labels.size(); iScatter++) {
        for (size_t iPoint = 0; iPoint < this->scatter_labels[iScatter].size(); iPoint++) {
            //                this->graphicLabels[iScatter][iPoint]->setPos(QPointF(this->scatter_data[iScatter][iPoint].first, this->scatter_data[iScatter][iPoint].second)); // this->chartView->chart()->mapToPosition(QPointF(this->scatter_data[iScatter][iPoint].first, this->scatter_data[iScatter][iPoint].second)));
            this->graphicLabels[iScatter][iPoint]->setPos(this->chartView->chart()->mapToPosition(QPointF(this->scatter_data[iScatter][iPoint].x, this->scatter_data[iScatter][iPoint].y)));
        }
    }
    if (!this->selectedPlotData.empty() || !this->selectedScatterData.empty()) {
        this->draw();
    }
    //    this->blockSignals(false);
    return this;
}

Plotter* Plotter::selectData(const Vector3& pos)
{
    if (!pos.isValid()) return this;

    float minDist = 0.05f;
    this->selectedScatterData.clear();
    this->selectedPlotData.clear();

    if (pos.isValid()) {
        for (size_t i = 0; i < this->plot_data.size(); i++) {
            for (size_t j = 0; j < this->plot_data[i].size(); j++) {
                if ((plot_data[i][j] - pos).norm2() < minDist*minDist)
                    this->selectedPlotData.push_back({i, j});
            }
        }
        for (size_t i = 0; i < this->scatter_data.size(); i++) {
            for (size_t j = 0; j < this->scatter_data[i].size(); j++) {
                if ((scatter_data[i][j] - pos).norm2() < minDist*minDist)
                    this->selectedScatterData.push_back({i, j});
            }
        }
    }

    if (!this->selectedPlotData.empty() || !this->selectedScatterData.empty()) {
        this->chartView->lockView();
    } else {
        this->chartView->unlockView();
    }

    if (!this->displayedImage.empty()) {
        Q_EMIT clickedOnImage(pos * displayedImage.getDimensions(), this->displayedImage(pos * displayedImage.getDimensions()));
    }
    return this;
}

Plotter* Plotter::displayInfoUnderMouse(const Vector3 &relativeMousePos)
{
    if (this->displayedImage.empty() || relativeMousePos.minComp() < 0.f || relativeMousePos.maxComp() > 1.f)
        return this;
    std::ostringstream oss;
    Vector3 size = displayedImage.getDimensions();
    Vector3 position = relativeMousePos * size;
    Vector3 value = this->displayedImage(position);
    oss << "Mouse pos: " << int(position.x) << ", " << int(position.y) << " -- Value : (" << value.x() << ", " << value.y() << ", " << value.z() << ") ";
    this->mouseInfoLabel->setText(QString::fromStdString(oss.str()));
    return this;
}
*/











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

PlotModel *PlotModel::addPlot(std::vector<Vector3> data, std::string name, QColor color)
{
    this->plot_data.push_back(data);
    this->plot_names.push_back(name);
    this->plot_colors.push_back(color);
    return this;
}

PlotModel *PlotModel::addScatter(std::vector<Vector3> data, std::string name, std::vector<std::string> labels, std::vector<QColor> colors)
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

PlotModel *PlotModel::addImage(GridV3 image, bool clamped, bool normalized, bool absolute, Vector3 minColors, Vector3 maxColors)
{
    this->displayedImage = image;
    //    image = image.flip(false, true);
    if (image.empty()) return this;
    if (clamped) {
        float min = std::numeric_limits<float>::max();
        float max = std::numeric_limits<float>::lowest();
        image.iterate([&](size_t i) {
            min = std::min(min, image[i].minComp());
            max = std::max(max, image[i].maxComp());
        });
        image.iterateParallel([&](size_t i) {
            for (int c = 0; c < 3; c++) {
                image[i][c] = std::clamp(image[i][c], minColors[c], maxColors[c]);
            }
        });
    }
    if (absolute) {
        image = image.abs();
    }
    if (normalized) {
        for (int c = 0; c < 3; c++) {
            float min = std::numeric_limits<float>::max();
            float max = std::numeric_limits<float>::lowest();
            image.iterate([&](size_t i) {
                min = std::min(image[i][c], min);
                max = std::max(image[i][c], max);
            });
            float d = max - min;
            if (d == 0) {
                image.iterateParallel([&](size_t i) {
                    image[i][c] = 0.f;
                });
            } else {
                image.iterateParallel([&](size_t i) {
                    image[i][c] = (image[i][c] - min) / d;
                });
            }
        }
        //        image.normalize();
    }
    Vector3 colorFilter = Vector3(1.f, 1.f, 1.f); //Vector3((displayR ? 1.f : 0.f), (displayG ? 1.f : 0.f), (displayB ? 1.f : 0.f));
    image.iterateParallel([&](size_t i) {
        image[i] *= colorFilter;
    });
    unsigned char* data = new unsigned char[image.size() * 4];

    for (size_t i = 0; i < image.size(); ++i) {
        data[int(4 * i + 2)] = (unsigned char)(std::clamp(image[i].x(), 0.f, 1.f) * 255);
        data[int(4 * i + 1)] = (unsigned char)(std::clamp(image[i].y(), 0.f, 1.f) * 255);
        data[int(4 * i + 0)] = (unsigned char)(std::clamp(image[i].z(), 0.f, 1.f) * 255);
        data[int(4 * i + 3)] = (unsigned char) 255;       // Alpha
    }

    if (this->backImage) {
        delete this->backImage;
    }
    this->backImage = new QImage(data, image.sizeX, image.sizeY, QImage::Format_ARGB32);
    //    *(this->backImage) = this->backImage->mirrored();
    return this;
}

PlotModel* PlotModel::reset()
{
    this->backImage = nullptr;
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
    return this;
}
