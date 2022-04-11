#ifndef CUSTOMINTERACTIVEOBJECT_H
#define CUSTOMINTERACTIVEOBJECT_H

#include <QObject>
#include <QWidget>

class CustomInteractiveObject : public QWidget
{
    Q_OBJECT
public:
    explicit CustomInteractiveObject(QWidget *parent = nullptr);

    virtual void hide() { this->visible = false; }
    virtual void show() { this->visible = true; }

    bool eventFilter(QObject* obj, QEvent* event);

    bool isHidden() { return !this->isVisible(); }
    bool isVisible() {
        if (!__initialized) this->setVisibility(this->visible);
        return this->visible;
    }
    bool setVisibility(bool visible) {
        this->visible = visible;
        if (visible) show();
        else hide();
        this->__initialized = true;
        return visible; }

Q_SIGNALS:
    void update();

protected:
    bool visible = true;
    bool __initialized = false;
    std::string id_name;
};

#endif // CUSTOMINTERACTIVEOBJECT_H
