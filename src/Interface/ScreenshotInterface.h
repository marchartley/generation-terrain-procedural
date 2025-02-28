#ifndef SCREENSHOTINTERFACE_H
#define SCREENSHOTINTERFACE_H

#include "ActionInterface.h"

class ScreenshotInterface : public ActionInterface
{
    Q_OBJECT
public:
    ScreenshotInterface(QWidget *parent = nullptr);

    void display(const Vector3& camPos = Vector3(false));

    void replay(nlohmann::json action);

    // void mouseMoveEvent(QMouseEvent* event);
    // void keyPressEvent(QKeyEvent* event);
    // void keyReleaseEvent(QKeyEvent* event);
    // void wheelEvent(QWheelEvent* event);
    // void mousePressEvent(QMouseEvent* event);

    QLayout* createGUI();

public Q_SLOTS:
    void takeScreenshots();
};

#endif // SCREENSHOTINTERFACE_H
