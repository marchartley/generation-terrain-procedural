#ifndef IMAGEVIEWER_H
#define IMAGEVIEWER_H

#include "GUIElements/AbstractPlotter.h"

class ImageViewer : public AbstractPlotter {
    Q_OBJECT
protected: // Singleton
    ImageViewer(const std::string& name, QWidget* parent = nullptr);
    ImageViewer(const std::string& name, ChartView* chartView, QWidget* parent = nullptr);

public:
    DECLARE_PLOTTER_GETTER(ImageViewer)

public Q_SLOTS:
    ImageViewer* updateUI();

    ImageViewer* setNormalizedModeImage(bool normalize);
    ImageViewer* setAbsoluteModeImage(bool absolute);
    ImageViewer* setFilteredValuesImage(bool filtered);

    ImageViewer* updateViewOptionsInterface();

    ImageViewer* displayInfoUnderMouse(const Vector3 &relativeMousePos);
};

#endif // IMAGEVIEWER_H
