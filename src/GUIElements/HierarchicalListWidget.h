#ifndef HIERARCHICALLISTWIDGET_H
#define HIERARCHICALLISTWIDGET_H

#include <QListWidget>
#include <QListWidgetItem>

class HierarchicalListWidgetItemBase;
template <class T>
class HierarchicalListWidgetItem;

enum HIERARCHY_TYPE {
    CHILD =0,
    SIBLING = 1,
    PARENT =2
};

class HierarchicalListWidget : public QListWidget
{
    Q_OBJECT
public:
    HierarchicalListWidget(QWidget* parent = nullptr);

    void setCurrentItem(int indexToSelect);
    void setCurrentItems(std::vector<int> indicesToSelect);

    template <class T>
    void removeItem(const T& toRemove)
    {
        for (int i = this->count() - 1; i >= 0; i--) {
            if (auto asTypedItem = dynamic_cast<HierarchicalListWidgetItem<T>*>(this->item(i))) {
                if (asTypedItem->itemValue == toRemove) {
                    this->takeItem(i);
                }
            }
        }
    }

Q_SIGNALS:
    void itemChangedHierarchy(int changedItemID, int relationItemID, HIERARCHY_TYPE newRelation, QDropEvent* event = nullptr);

public Q_SLOTS:
    void dragEnterEvent(QDragEnterEvent* event);
    void dragLeaveEvent(QDragLeaveEvent* event);
    void dropEvent(QDropEvent* event);

public:
    HierarchicalListWidgetItemBase* movingItem = nullptr;
};

class HierarchicalListWidgetItemBase : public QListWidgetItem
{
    //    Q_OBJECT
public:
    HierarchicalListWidgetItemBase(const std::string& internal_text = "", QListWidget* parent = nullptr);
    HierarchicalListWidgetItemBase(const std::string& internal_text, int ID, int level = 0, QListWidget* parent = nullptr);

    void setLevel(int newLevel);

    int ID;
    int level;
    std::string internalText;
};

template <class T = int>
class HierarchicalListWidgetItem : public HierarchicalListWidgetItemBase
{
//    Q_OBJECT
public:
    HierarchicalListWidgetItem(const std::string& internal_text, T itemValue = T(), int ID = 0, int level = 0, QListWidget* parent = nullptr)
        : HierarchicalListWidgetItemBase(internal_text, ID, level, parent), itemValue(itemValue)
    {}

    T getItem() const { return this->itemValue; }
    void setItem(T item) { this->itemValue = item; }

    T itemValue;
};

#endif // HIERARCHICALLISTWIDGET_H
