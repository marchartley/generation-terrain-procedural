#ifndef CUSTOMINTERACTIVEOBJECT_H
#define CUSTOMINTERACTIVEOBJECT_H

#include <QObject>

class CustomInteractiveObject : public QObject
{
    Q_OBJECT
public:
    explicit CustomInteractiveObject();

    virtual void hide() { this->visible = false; }
    virtual void show() { this->visible = true; }
    bool isHidden() { return !this->isVisible(); }
    bool isVisible() { return this->visible; }
    bool setVisibility(bool visible) {
        this->visible = visible;
        if (visible) show();
        else hide();
        return visible; }

Q_SIGNALS:
    void update();

protected:
    bool visible = true;
};

#endif // CUSTOMINTERACTIVEOBJECT_H
