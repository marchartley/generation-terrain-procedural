#include "StickyFrame.h"

#include <iostream>
#include <QHBoxLayout>
#include <QPushButton>


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

StickyFrame::StickyFrame(QWidget *parent, float x, float y, float w, float h, bool useAbsolutePosition)
    : StickyFrame(parent)
{
    parent->installEventFilter(this);
    this->x = x;
    this->y = y;
    this->w = w;
    this->h = h;
    this->useAbsolutePosition = useAbsolutePosition;
}

void StickyFrame::resizeEvent(QResizeEvent *e)
{
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
    content->setLayout(contentLayout);
    return this;
}

void StickyFrame::resizeWithParent()
{
    if (this->useAbsolutePosition) {
        this->setGeometry(x, y, w, h);
    } else {
        this->setGeometry(
                    this->parentWidget()->width() * x,
                    this->parentWidget()->height() * y,
                    this->parentWidget()->width() * w,
                    this->parentWidget()->height() * h
                    );
    }
}
