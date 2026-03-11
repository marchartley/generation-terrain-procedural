#ifndef GANUIPAINTER_H
#define GANUIPAINTER_H

#include "GUIElements/AbstractPlotter.h"

class GanUIPainter : public AbstractPlotter {
    Q_OBJECT
public: // private: // Singleton
    GanUIPainter(const std::string& name, QWidget* parent = nullptr);

public:
    DECLARE_PLOTTER_GETTER(GanUIPainter)

    void drawStroke(const Vector3& pStart, const Vector3& pEnd);

    static Vector3 getColorFromIndex(float index);

    Vector3 currentColor = Vector3(1.f, .5f, .5f);

    float brushSize = 50.f;
    float colorIndex = .2f;
    float sharpness = 0.01f;

    GanUIPainter& updateUI();
public Q_SLOTS:

    // Q_SIGNALS:
    // void clickedOnImage(const Vector3& pos, Vector3 value);
    // void movedOnImage(const Vector3& pos, const Vector3& previousPos, QMouseEvent* event);
};

#endif // GANUIPAINTER_H
