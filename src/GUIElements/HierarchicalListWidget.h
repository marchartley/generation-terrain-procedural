#ifndef HIERARCHICALLISTWIDGET_H
#define HIERARCHICALLISTWIDGET_H

#include <QListWidget>
#include <QListWidgetItem>

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

Q_SIGNALS:
    void itemChangedHierarchy(int changedItemID, int relationItemID, HIERARCHY_TYPE newRelation, QDropEvent* event = nullptr);

public Q_SLOTS:
    void dragEnterEvent(QDragEnterEvent* event);
    void dragLeaveEvent(QDragLeaveEvent* event);
    void dropEvent(QDropEvent* event);

public:
    HierarchicalListWidgetItem* movingItem = nullptr;
};

class HierarchicalListWidgetItem : public QListWidgetItem
{
//    Q_OBJECT
public:
    HierarchicalListWidgetItem(std::string internal_text = "", QListWidget* parent = nullptr);
    HierarchicalListWidgetItem(std::string internal_text, int ID, int level = 0, QListWidget* parent = nullptr);

    void setLevel(int newLevel);


    int ID;
    int level;
    std::string internalText;
};

#endif // HIERARCHICALLISTWIDGET_H
