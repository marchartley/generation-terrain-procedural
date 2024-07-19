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
}

Vector3 ChartView::getRelativeMousePositionInImage(const Vector3 &pos)
{
    // Get the coordinate in the plotted area...
    QPointF qMousePos(pos.x, pos.y);
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
Plotter* Plotter::instance = nullptr;

Plotter::Plotter(QWidget *parent) : Plotter(new ChartView(new Chart()), parent)
{
}

Plotter::Plotter(ChartView *chartView, QWidget *parent) : QDialog(parent), chartView(chartView)
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

Plotter *Plotter::getInstance()
{
    if (Plotter::instance == nullptr) {
//        std::cerr << "Plotter has not been initialized with function Plotter::init()" << std::endl;
        Plotter::init();
    }
    return Plotter::instance;
}

Plotter *Plotter::init(ChartView *chartView, QWidget *parent)
{
    Plotter::instance = new Plotter(chartView, parent);
    return Plotter::getInstance();
}

Plotter* Plotter::addPlot(std::vector<float> data, std::string name, QColor color)
{
    std::vector<Vector3> _data;
    for (unsigned int i = 0; i < data.size(); i++) {
        _data.push_back(Vector3(i, data[i]));
    }
    this->addPlot(_data, name, color);
    return Plotter::getInstance();
}

Plotter* Plotter::addPlot(std::vector<Vector3> data, std::string name, QColor color)
{
    this->plot_data.push_back(data);
    this->plot_names.push_back(name);
    this->plot_colors.push_back(color);
    return Plotter::getInstance();
}

Plotter* Plotter::addScatter(std::vector<float> data, std::string name, std::vector<std::string> labels, std::vector<QColor> colors)
{
    std::vector<Vector3> _data;
    for (unsigned int i = 0; i < data.size(); i++) {
        _data.push_back(Vector3(i, data[i]));
    }
    this->addScatter(_data, name, labels, colors);
    return Plotter::getInstance();
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
    return Plotter::getInstance();
}

Plotter* Plotter::addImage(GridV3 image)
{
    this->displayedImage = image;
//    image = image.flip(false, true);
    if (image.empty()) return Plotter::getInstance();
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
            image.iterateParallel([&](size_t i) {
                image[i][c] = (image[i][c] - min) / d;
            });
        }
//        image.normalize();
    }
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
    return Plotter::getInstance();
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
    imgSize.z = 1;
    GridV3 img(imgSize);
    Vector3 reducedSize = (field.getDimensions() * reductionFactor).roundedUp(); //(30, 30, 1);
    reducedSize.z = 1;
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
    imgSize.z = 1;
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
    return Plotter::getInstance();
}

Plotter* Plotter::show()
{
    this->draw();
    QDialog::show();
    return Plotter::getInstance();
}

Plotter* Plotter::updateUI()
{
    blockSignals(true);
    /*for (auto& child : this->children()) {
        if (auto asUIElement = dynamic_cast<UIElement*>(child)) {
            asUIElement->update();
        }
    }*/
    this->interfaceButtons->update();
    /*for (auto& child : this->interfaceButtons->children()) {
        if (auto asUIElement = dynamic_cast<UIElement*>(child)) {
            asUIElement->update();
        }
    }*/
    blockSignals(false);
    return Plotter::getInstance();
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
    return Plotter::getInstance();
}

Plotter* Plotter::copyToClipboard()
{
    QPixmap p = this->chartView->grab();
    QApplication::clipboard()->setPixmap(p, QClipboard::Clipboard);
    return Plotter::getInstance();
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

    return Plotter::getInstance();
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
    return Plotter::getInstance();
}

Plotter* Plotter::selectData(const Vector3& pos)
{
    if (!pos.isValid()) return Plotter::getInstance();

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
    return Plotter::getInstance();
}

Plotter* Plotter::displayInfoUnderMouse(const Vector3 &relativeMousePos)
{
    if (this->displayedImage.empty() || relativeMousePos.minComp() < 0.f || relativeMousePos.maxComp() > 1.f)
        return Plotter::getInstance();
    std::ostringstream oss;
    Vector3 size = displayedImage.getDimensions();
    Vector3 position = relativeMousePos * size;
    Vector3 value = this->displayedImage(position);
    oss << "Mouse pos: " << int(position.x) << ", " << int(position.y) << " -- Value : (" << value.x << ", " << value.y << ", " << value.z << ") ";
    this->mouseInfoLabel->setText(QString::fromStdString(oss.str()));
    return Plotter::getInstance();
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

