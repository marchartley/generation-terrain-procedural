#ifndef GANUIPAINTER_H
#define GANUIPAINTER_H

#include "Graphics/AbstractPlotter.h"

class GanUIPainter : public AbstractPlotter {
    Q_OBJECT
private: // Singleton
    GanUIPainter(std::string name, QWidget* parent = nullptr);
    GanUIPainter(std::string name, ChartView* chartView, QWidget* parent = nullptr);

public:
    static GanUIPainter* getInstance(std::string name = "");
    static GanUIPainter* get(std::string name = "") { return GanUIPainter::getInstance(toLower(name)); }
    static GanUIPainter* init(std::string name, ChartView* chartView = nullptr, QWidget* parent = nullptr);

    void drawStroke(const Vector3& pStart, const Vector3& pEnd);

    static Vector3 getColorFromIndex(float index);

    Vector3 currentColor = Vector3(1.f, .5f, .5f);

    float brushSize = 50.f;
    float colorIndex = .2f;
    float sharpness = 0.01f;

public Q_SLOTS:
    GanUIPainter* updateUI();

    // Q_SIGNALS:
    // void clickedOnImage(const Vector3& pos, Vector3 value);
    // void movedOnImage(const Vector3& pos, const Vector3& previousPos, QMouseEvent* event);
};

#endif // GANUIPAINTER_H
