#include "StickyFrame.h"

#include <iostream>
#include <QHBoxLayout>
#include <QPushButton>
#include <QTimer>


StickyFrame::StickyFrame(QWidget* parent) : QWidget(parent)
{
    QHBoxLayout* layout = new QHBoxLayout(this);
    content = new QFrame();
    content->setAutoFillBackground(true);
    QPalette pal = QPalette();

    // set black background
    // Qt::black / "#000000" / "black"
    pal.setColor(QPalette::Window, Qt::white);
    content->setPalette(pal);
    /*this->setStyleSheet(QString("QFrame{background-color: rgba(1.0, 1.0, 1.0, 0.5); "
                                   "color : rgba(1.0, 1.0, 1.0, 1.0);"
                                   "}"
                                   "QFrame * {color : white; }"));*/
    content->setLayout(new QVBoxLayout());
    layout->addWidget(content);
    this->setLayout(layout);
}

StickyFrame::StickyFrame(QWidget *parent, float x0, float y0, float x1, float y1, bool useAbsolutePosition)
    : StickyFrame(parent)
{
    parent->installEventFilter(this);
    this->x0 = x0;
    this->y0 = y0;
    this->x1 = x1;
    this->y1 = y1;
    this->useAbsolutePosition = useAbsolutePosition;
}

void StickyFrame::resizeEvent(QResizeEvent *e)
{
    QWidget::resizeEvent(e);
}

bool StickyFrame::eventFilter(QObject *obj, QEvent *event)
{
    if (event->type() == QEvent::Resize) {
        this->resizeWithParent();
    }
    return false; //QObject::eventFilter(obj, event);
}

void StickyFrame::clearContent()
{
    this->content->setUpdatesEnabled(false);
    qDeleteAll(content->findChildren<QWidget*>("", Qt::FindDirectChildrenOnly));
    this->content->setUpdatesEnabled(true);
}

StickyFrame* StickyFrame::setContent(QLayout* contentLayout)
{
    qDeleteAll(content->findChildren<QWidget *>(QString(), Qt::FindDirectChildrenOnly));
    delete content->layout();
    content->setSizePolicy(QSizePolicy::Policy::Expanding, QSizePolicy::Policy::Expanding);
    content->setLayout(contentLayout);
    content->layout()->activate();
    this->layout()->activate();
    resizeWithParent();

    // Worst solution maybe... but fuck it, it works!
    QTimer::singleShot(10, this, &StickyFrame::resizeWithParent);
    return this;
}

void StickyFrame::resizeWithParent()
{
    float x0 = this->x0;
    float y0 = this->y0;
    float x1 = this->x1;
    float y1 = this->y1;

    this->adjustSize();
    float estimatedWidth = this->sizeHint().width();
    float estimatedHeight = this->sizeHint().height();

    if (this->useAbsolutePosition) {

    } else {
        if (x0 == -1) {
            if (x1 == -1) {
                // Just center the widget
                x0 = this->parentWidget()->width() - estimatedWidth/2;
                x1 = this->parentWidget()->width() + estimatedWidth/2;
            } else {
                x1 = this->parentWidget()->width() * x1;
                x0 = std::max(x1 - estimatedWidth, 0.f);
            }
        } else {
            x0 = this->parentWidget()->width() * x0;
            if (x1 == -1) {
                x1 = x0 + estimatedWidth;
            } else {
                x1 = this->parentWidget()->width() * x1;
            }
        }

        if (y0 == -1) {
            if (y1 == -1) {
                // Just center the widget
                y0 = this->parentWidget()->height() - estimatedHeight/2;
                y1 = this->parentWidget()->height() + estimatedHeight/2;
            } else {
                y1 = this->parentWidget()->height() * y1;
                y0 = std::max(y1 - estimatedHeight, 0.f);
            }
        } else {
            y0 = this->parentWidget()->height() * y0;
            if (y1 == -1) {
                y1 = y0 + estimatedHeight;
            } else {
                y1 = this->parentWidget()->height() * y1;
            }
        }
    }
    this->setGeometry(x0, y0, x1 - x0, y1 - y0);
}
