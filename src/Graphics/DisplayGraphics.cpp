#include "DisplayGraphics.h"
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

bool ChartView::viewportEvent(QEvent *event)
{
    return QChartView::viewportEvent(event);
}

void ChartView::mousePressEvent(QMouseEvent *event)
{
    this->previousMousePos = event->pos();

    QPointF local = this->chart()->mapToValue(this->previousMousePos);
    Vector3 vecLocal = Vector3(local.x(), local.y());

    Q_EMIT this->clickedOnValue(vecLocal);

    return QChartView::mousePressEvent(event);
}
void ChartView::mouseMoveEvent(QMouseEvent *event)
{
    if (event->buttons().testFlag(Qt::LeftButton)) {
        QPoint delta = event->pos() - this->previousMousePos;
        this->previousMousePos = event->pos();
        if (!this->locked)
            this->chart()->scroll(-delta.x(), delta.y());
        Q_EMIT this->updated();
    }
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

Plotter::Plotter(QWidget *parent) : Plotter(new ChartView(new Chart()), parent)
{
}

Plotter::Plotter(ChartView *chartView, QWidget *parent) : QDialog(parent), chartView(chartView)
{
    this->setLayout(new QHBoxLayout());
    this->layout()->addWidget(this->chartView);
    chartView->setRenderHint(QPainter::Antialiasing);
    this->chartView->chart()->legend()->setMarkerShape(QLegend::MarkerShapeFromSeries);
    this->resize(800, 600);

//    this->m_axisX = new QValueAxis(this->chartView->chart());
//    this->m_axisY = new QValueAxis(this->chartView->chart());
//    this->chartView->chart()->addAxis(m_axisX,Qt::AlignBottom);
//    this->chartView->chart()->addAxis(m_axisY,Qt::AlignLeft);

    this->saveButton = new QPushButton("Save");
    this->layout()->addWidget(this->saveButton);

    QObject::connect(this->saveButton, &QPushButton::pressed, [&]() {
        QString q_filename = QFileDialog::getSaveFileName(this, QString("Save plot"));
        if (!q_filename.isEmpty())
            saveFig(q_filename.toStdString());
    });
    QObject::connect(this->chartView, &ChartView::clickedOnValue, this, &Plotter::selectData);
    QObject::connect(this->chartView->chart(), &QChart::geometryChanged, this, &Plotter::updateLabelsPositions);
    QObject::connect(this->chartView, &ChartView::updated, this, &Plotter::updateLabelsPositions);
}

void Plotter::addPlot(std::vector<float> data, std::string name, QColor color)
{
    std::vector<Vector3> _data;
    for (unsigned int i = 0; i < data.size(); i++) {
        _data.push_back(Vector3(i, data[i]));
    }
    this->addPlot(_data, name, color);
}

void Plotter::addPlot(std::vector<Vector3> data, std::string name, QColor color)
{
    this->plot_data.push_back(data);
    this->plot_names.push_back(name);
    this->plot_colors.push_back(color);
}

void Plotter::addScatter(std::vector<float> data, std::string name, std::vector<std::string> labels, std::vector<QColor> colors)
{
    std::vector<Vector3> _data;
    for (unsigned int i = 0; i < data.size(); i++) {
        _data.push_back(Vector3(i, data[i]));
    }
    this->addScatter(_data, name, labels, colors);
}

void Plotter::addScatter(std::vector<Vector3> data, std::string name, std::vector<std::string> labels, std::vector<QColor> colors)
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
}

void Plotter::draw()
{
//    QTransform prevState = this->chartView->transform();
    auto prevState = this->chartView->chart()->transformations();
    this->chartView->chart()->removeAllSeries();
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

        /*for (size_t j = 0; j < this->scatter_data[i].size(); j++) {
            if (PlotColorToQColor.count(this->scatter_colors[i][j])) {
                QGraphicsItem *it = this->chartView->itemAt(this->chartView->mapFromScene(this->chartView->chart()->mapToPosition(series->at(j))));
                if(QGraphicsEllipseItem *ellipse = qgraphicsitem_cast<QGraphicsEllipseItem*>(it)){
                    QColor color = QColor::fromRgb(QRandomGenerator::global()->generate());;
                    ellipse->setBrush(color);
                }
            }
        }*/
        if (series->name().isEmpty()) {
            this->chartView->chart()->legend()->markers(series)[0]->setVisible(false);
        }
    }
    while (!this->chartView->chart()->axes().empty()) {
        this->chartView->chart()->removeAxis(this->chartView->chart()->axes().front());
    }
    this->chartView->chart()->createDefaultAxes();
    if (!title.empty())
        this->chartView->chart()->setTitle(QString::fromStdString(title));

    for (auto& labels : this->graphicLabels)
        for (auto& lab : labels)
            delete lab;
    this->graphicLabels.clear();
    for (size_t iScatter = 0; iScatter < this->scatter_labels.size(); iScatter++) {
        this->graphicLabels.push_back(std::vector<QGraphicsTextItem*>());
        for (size_t iPoint = 0; iPoint < this->scatter_labels[iScatter].size(); iPoint++) {
            QGraphicsTextItem *itm = new QGraphicsTextItem(QString::fromStdString(this->scatter_labels[iScatter][iPoint]), this->chartView->chart());
            this->graphicLabels[iScatter].push_back(itm);
        }
    }
//    this->chartView->setTransform(prevState);
//    this->chartView->chart()->setTransformations(prevState);
    this->chartView->chart()->zoomOut();
    this->chartView->update();
}

int Plotter::exec()
{
    this->draw();

    return QDialog::exec();

}

void Plotter::saveFig(std::string filename)
{
    QPixmap p = this->chartView->grab();
    p.save(QString::fromStdString(filename), "PNG");
}

void Plotter::updateLabelsPositions()
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
}

void Plotter::selectData(Vector3 pos)
{
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
