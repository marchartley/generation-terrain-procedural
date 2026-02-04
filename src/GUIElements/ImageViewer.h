#ifndef IMAGEVIEWER_H
#define IMAGEVIEWER_H

#include "GUIElements/AbstractPlotter.h"

class ImageViewer : public AbstractPlotter {
    Q_OBJECT
protected: // Singleton
    ImageViewer(const std::string& name, QWidget* parent = nullptr);
    ImageViewer(const std::string& name, ChartView* chartView, QWidget* parent = nullptr);

public:
    static ImageViewer* getInstance(std::string name = "");
    static ImageViewer* get(std::string name = "") { return ImageViewer::getInstance(toLower(name)); }
    static ImageViewer* init(const std::string& name, ChartView* chartView = nullptr, QWidget* parent = nullptr);

    // RangeSliderElement* rangeValuesWidget;

public Q_SLOTS:
    ImageViewer* updateUI();

    ImageViewer* setNormalizedModeImage(bool normalize);
    ImageViewer* setAbsoluteModeImage(bool absolute);
    ImageViewer* setFilteredValuesImage(bool filtered);

    ImageViewer* updateViewOptionsInterface();

// Q_SIGNALS:
    // void clickedOnImage(const Vector3& pos, Vector3 value);
    // void movedOnImage(const Vector3& pos, const Vector3& previousPos, QMouseEvent* event);
};

#endif // IMAGEVIEWER_H
