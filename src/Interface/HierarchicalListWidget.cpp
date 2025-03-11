#include "HierarchicalListWidget.h"
#include <iostream>
#include <QDropEvent>

#include "Utils/Utils.h"

HierarchicalListWidget::HierarchicalListWidget(QWidget *parent)
    : QListWidget(parent)
{

    this->setDragDropMode(QAbstractItemView::DragDrop);
    this->setDefaultDropAction(Qt::MoveAction);
    this->setSelectionMode(QAbstractItemView::SingleSelection);
}

void HierarchicalListWidget::setCurrentItem(int indexToSelect)
{
    this->clearSelection();
    for (auto& _child : this->findItems("*", Qt::MatchWildcard)) {
        HierarchicalListWidgetItem* item = dynamic_cast<HierarchicalListWidgetItem*>(_child);
        if (item != nullptr) {
            if (item->ID == indexToSelect) {
                item->setSelected(true);
                return;
            }
        }
    }
}

void HierarchicalListWidget::setCurrentItems(std::vector<int> indicesToSelect)
{
    QObject::blockSignals(true);
    this->clearSelection();
    for (auto& _child : this->findItems("*", Qt::MatchWildcard)) {
        HierarchicalListWidgetItem* item = dynamic_cast<HierarchicalListWidgetItem*>(_child);
        if (item != nullptr) {
            if (isIn(item->ID, indicesToSelect)) {
                item->setSelected(true);
            }
        }
    }
    QObject::blockSignals(false);
}

void HierarchicalListWidget::dragEnterEvent(QDragEnterEvent *event)
{

    if (this->selectedItems().empty()) return;

    this->movingItem = dynamic_cast<HierarchicalListWidgetItem*>(this->selectedItems()[0]);
    if (movingItem != nullptr) {
        std::cout << "Dragging item #" << this->movingItem->ID << std::endl;
    }

    QListWidget::dragEnterEvent(event);
}

void HierarchicalListWidget::dragLeaveEvent(QDragLeaveEvent *event)
{
    QListWidget::dragLeaveEvent(event);
}

void HierarchicalListWidget::dropEvent(QDropEvent *event)
{
    /*
    if (movingItem != nullptr) {

        std::cout << (this->dropIndicatorPosition() == QAbstractItemView::AboveItem ? "above" : (this->dropIndicatorPosition() == QAbstractItemView::BelowItem ? "below" : (this->dropIndicatorPosition() == QAbstractItemView::OnItem ? "exact" : ""))) << std::endl;
    }
    this->movingItem = nullptr;
    */
    bool droppedAbove = this->dropIndicatorPosition() == QAbstractItemView::AboveItem;
    bool droppedBelow = this->dropIndicatorPosition() == QAbstractItemView::BelowItem;

    if (movingItem != nullptr) {
        int droppedRow = this->indexAt(event->pos()).row();
        std::cout << "Dropped on row " << droppedRow << std::endl;
        HierarchicalListWidgetItem* linkedItem = dynamic_cast<HierarchicalListWidgetItem*>(this->item(droppedRow));
        if (droppedAbove) {
            QListWidget::dropEvent(event);
            Q_EMIT this->itemChangedHierarchy(movingItem->ID, linkedItem->ID, SIBLING, event);
//            int previousLevel = movingItem->level;
//            movingItem->setLevel(linkedItem->level);
            /*for (int i = this->currentRow(); i < this->model()->rowCount(); i++) {
                auto child = dynamic_cast<HierarchicalListWidgetItem*>(this->item(i));
                if (child->level > previousLevel) {
                    child->
                } else {
                    break;
                }
            }*/
        } else if (droppedBelow) {
            QListWidget::dropEvent(event);
            Q_EMIT this->itemChangedHierarchy(movingItem->ID, linkedItem->ID, CHILD, event);
//            movingItem->setLevel(linkedItem->level + 1);
        } else {
            std::cout << "Wasn't meant to be droppe here...";
        }
        this->movingItem = nullptr;
    }
}

HierarchicalListWidgetItem::HierarchicalListWidgetItem(std::string internal_text, QListWidget *parent)
    : HierarchicalListWidgetItem(internal_text, -1, 0, parent)
{

}

HierarchicalListWidgetItem::HierarchicalListWidgetItem(std::string internal_text, int ID, int level, QListWidget *parent)
    : QListWidgetItem(parent), ID(ID), level(level), internalText(internal_text)
{
    this->setLevel(level);
}

void HierarchicalListWidgetItem::setLevel(int newLevel)
{
    this->level = std::max(newLevel, 0);
    this->setText(QString(level * 2, ' ') + QString::fromStdString(this->internalText));
}
