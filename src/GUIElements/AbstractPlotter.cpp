#include "AbstractPlotter.h"


#include <iostream>







std::map<std::string, AbstractPlotter*> AbstractPlotter::instances = std::map<std::string, AbstractPlotter*>();
std::string AbstractPlotter::defaultName = "default";

AbstractPlotter::AbstractPlotter(const std::string &name, QWidget *parent)
    : AbstractPlotter(name, name, parent)
{
}

AbstractPlotter::AbstractPlotter(const std::string& name, const std::string &title, QWidget *parent) : QDialog(parent), name(name)
{
    this->chartView = new ChartView();

    this->dataModel = std::make_shared<PlotModel>();

    auto layout = new InterfaceUI(InterfaceUI::VERTICAL);
    auto mainLayout = new InterfaceUI(InterfaceUI::HORIZONTAL);

    mainInterface = new InterfaceUI(InterfaceUI::VERTICAL, "Main");
    mainInterface->add(new UIElement(this->chartView));
    toolsInterface = new InterfaceUI(InterfaceUI::VERTICAL, "Tools");
    viewOptionsInterface = new InterfaceUI(InterfaceUI::VERTICAL, "View options");
    saveCopyInterface = new InterfaceUI(InterfaceUI::VERTICAL, "Save/Copy");
    infosInterface = new InterfaceUI(InterfaceUI::HORIZONTAL, "Infos");

    this->chartView->setRenderHint(QPainter::Antialiasing);
    this->chartView->chart()->legend()->setMarkerShape(QLegend::MarkerShapeFromSeries);
    viewAndCopyInterface = new InterfaceUI(InterfaceUI::VERTICAL);
    viewAndCopyInterface->add({viewOptionsInterface, saveCopyInterface});
    //    this->chartView->setMaximumSize(10000, 10000);
    //    this->chartView->chart()->setMaximumSize(10000, 10000);
    //    this->chartView->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    //    this->chartView->chart()->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    //    this->mouseInfoLabel->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Maximum);


    layout->add({mainLayout, infosInterface});

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

    saveCopyInterface->add({saveButton, copyToClipboardButton});

    this->setLayout(layout->get()->layout());

    this->resize(1100, 600);
    this->updateGeometry();

    this->setWindowTitle(QString::fromStdString(toCapitalize(title)));

    // QObject::connect(this->chartView, &ChartView::clickedOnValue, this, &AbstractPlotter::selectData);
    QObject::connect(this->chartView, &ChartView::mouseMoved, this, [&](const Vector3& pos, const Vector3& prevPos, QMouseEvent* e){
        this->displayInfoUnderMouse(pos);
        if (this->hasImage()) {
            Vector3 scale = dataModel->getImage().getDimensions();
            this->emitMouseMoved(pos * scale, prevPos * scale, e);
            emitMouseMoved(pos * scale, prevPos * scale, e);
        }
        else if (this->hasVectorField()) {
            Vector3 scale = dataModel->vectorData.field.getDimensions();
            this->emitMouseMoved(pos * scale, prevPos * scale, e);
            emitMouseMoved(pos * scale, prevPos * scale, e);
        }
    });
    // QObject::connect(this->chartView->chart(), &QChart::geometryChanged, this, &AbstractPlotter::draw);
    QObject::connect(this->chartView, &ChartView::clickedOnValue, this, [&](const Vector3& pos, bool leftClick, bool rightClick) {
        if (pos.isValid()) {
            Vector3 imageValue = (this->hasImage() ? this->dataModel->getImage().at(pos * this->dataModel->getImage().getDimensions()) : Vector3::invalid);
            emitMousePressed(pos, imageValue, leftClick, rightClick);
            // Q_EMIT this->clickedOnImage(pos, imageValue, leftClick, rightClick);
        } else {
            emitMouseReleased(pos.validated());
        }
    });
}

AbstractPlotter& AbstractPlotter::addPlot(const std::vector<float>& data, std::string name, QColor color)
{
    std::vector<Vector3> _data;
    for (unsigned int i = 0; i < data.size(); i++) {
        _data.push_back(Vector3(i, data[i]));
    }
    return this->addPlot(_data, name, color);
}

AbstractPlotter& AbstractPlotter::addPlot(const std::vector<Vector3>& data, std::string name, QColor color)
{
    this->dataModel->addPlot(data, name, color);
    return *this;
}

AbstractPlotter& AbstractPlotter::addPlot(const BSpline &data, std::string name, QColor color)
{
    return this->addPlot(data.getPoints(), name, color);
}

AbstractPlotter& AbstractPlotter::addScatter(const std::vector<float>& data, std::string name, std::vector<std::string> labels, std::vector<QColor> colors)
{
    std::vector<Vector3> _data;
    for (unsigned int i = 0; i < data.size(); i++) {
        _data.push_back(Vector3(i, data[i]));
    }
    return this->addScatter(_data, name, labels, colors);
}

AbstractPlotter &AbstractPlotter::addScatter(const std::vector<float>& dataX, const std::vector<float>& dataY, const std::string& name, const std::vector<std::string>& labels, const std::vector<QColor>& colors)
{
    std::vector<Vector3> _data;
    for (size_t i = 0; i < dataX.size(); i++) {
        _data.push_back(Vector3(dataX[i], dataY[i]));
    }
    return this->addScatter(_data, name, labels, colors);
}

AbstractPlotter& AbstractPlotter::addScatter(const std::vector<Vector3>& data, std::string name, std::vector<std::string> labels, std::vector<QColor> colors)
{
    this->dataModel->addScatter(data, name, labels, colors);
    return *this;
}

AbstractPlotter& AbstractPlotter::addImage(const GridV3& image)
{
    this->dataModel->addImage(image);
    return *this;
}

AbstractPlotter& AbstractPlotter::addImage(const GridF &image)
{
    this->dataModel->addImage(image);
    return *this;
}

AbstractPlotter& AbstractPlotter::addImage(const Matrix3<double> &image)
{
    return this->addImage((GridF)image);
}

AbstractPlotter& AbstractPlotter::addImage(const GridI &image)
{
    return this->addImage((GridF)image);
}

AbstractPlotter& AbstractPlotter::addVectorField(const GridV3 &field)
{
    this->dataModel->addVectorField(field);
    return *this;
}

AbstractPlotter& AbstractPlotter::setOverlay(const GridV3 &colors, const GridF &alpha, const std::string &overlayName, int layer)
{
    this->chartView->setOverlay(colors, overlayName, alpha, layer);
    return *this;
}

AbstractPlotter& AbstractPlotter::setOverlay(const std::pair<GridV3, GridF> &colorAndAlpha, const std::string &overlayName, int layer)
{
    return this->setOverlay(colorAndAlpha.first, colorAndAlpha.second, overlayName, layer);
}

AbstractPlotter& AbstractPlotter::setOverlay(const GridF& colors, const GridF &alpha, const std::string &overlayName, int layer)
{
    return this->setOverlay(this->dataModel->imageData.prepareImageForDisplay(Image(colors)), alpha, overlayName, layer);
}

AbstractPlotter& AbstractPlotter::setOverlay(const std::pair<GridF, GridF> &colorAndAlpha, const std::string &overlayName, int layer)
{
    return this->setOverlay(colorAndAlpha.first, colorAndAlpha.second, overlayName, layer);
}

AbstractPlotter& AbstractPlotter::showOverlay(const std::string &overlayName)
{
    this->chartView->overlayDisplayed[overlayName] = true;
    return *this;
}
AbstractPlotter& AbstractPlotter::hideOverlay(const std::string &overlayName)
{
    this->chartView->overlayDisplayed[overlayName] = false;
    return *this;
}

void AbstractPlotter::draw()
{
    this->chartView->setPlotModel(this->dataModel);
}

void AbstractPlotter::show()
{
    this->updateUI();
    this->draw();
    QDialog::show();
}

AbstractPlotter& AbstractPlotter::updateUI(bool forceUpdate)
{
    if (forceUpdate || !this->isVisible()) {
        blockSignals(true);

        this->updateToolsInterface();
        this->updateViewOptionsInterface();
        this->updateSaveCopyInterface();
        this->updateInfosInterface();
        blockSignals(false);
    }
    return *this;
}

AbstractPlotter& AbstractPlotter::updateToolsInterface()
{
    this->toolsInterface->update();
    return *this;
}

AbstractPlotter& AbstractPlotter::updateViewOptionsInterface()
{
    this->viewOptionsInterface->update();
    return *this;
}

AbstractPlotter& AbstractPlotter::updateSaveCopyInterface()
{
    this->saveCopyInterface->update();
    return *this;
}

AbstractPlotter& AbstractPlotter::updateInfosInterface()
{
    this->infosInterface->update();
    return *this;
}

int AbstractPlotter::exec()
{
    this->show();

    return QDialog::exec();

}

AbstractPlotter& AbstractPlotter::saveFig(const std::string& filename)
{
    QPixmap p = this->chartView->grab();
    if (this->hasImage())
        p = QPixmap::fromImage(this->dataModel->imageData.computeDisplayedImage());
    p.save(QString::fromStdString(filename), "PNG");
    std::cout << "Image " << filename << " saved." << std::endl;
    return *this;
}

AbstractPlotter& AbstractPlotter::copyToClipboard()
{
    QPixmap p = this->chartView->grab();
    QApplication::clipboard()->setPixmap(p, QClipboard::Clipboard);
    return *this;
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

void AbstractPlotter::hideEvent(QHideEvent *event)
{
    QDialog::hideEvent(event);
    // Flag the Chartview updater
}

QTimer *AbstractPlotter::animate(std::function<bool ()> callback, int interval_ms)
{
    QTimer* t = new QTimer(this);
    t->setInterval(interval_ms);
    t->setSingleShot(true);
    QObject::connect(t, &QTimer::timeout, [t, callback, interval_ms]() {
        if (callback())
            t->start(interval_ms);
    });
    t->start();
    return t;
}

/*
AbstractPlotter& AbstractPlotter::reset()
{
    this->dataModel->reset();
    this->chartView->setPlotModel(this->dataModel);
    return *this;
}
*/

void AbstractPlotter::displayInfoUnderMouse(const Vector3 &relativeMousePos)
{
    if (!this->hasImage() || relativeMousePos.minComp() < 0.f || relativeMousePos.maxComp() > 1.f)
        return;
    std::ostringstream oss;
    Vector3 size = this->dataModel->getImage().getDimensions();
    Vector3 position = relativeMousePos * size;
    Vector3 value = this->dataModel->getImage()(position);
    oss << "Mouse pos: " << int(position.x()) << ", " << int(position.y()) << " -- Value : (" << value.x() << ", " << value.y() << ", " << value.z() << ") [norm: " << value.norm() << "]";
    this->mouseInfoLabel->setText(QString::fromStdString(oss.str()));
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


