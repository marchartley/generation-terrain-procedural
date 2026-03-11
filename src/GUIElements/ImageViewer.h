#ifndef IMAGEVIEWER_H
#define IMAGEVIEWER_H

#include "GUIElements/AbstractPlotter.h"

class ImageViewer : public AbstractPlotter {
    Q_OBJECT
public: // protected: // Singleton
    ImageViewer(const std::string& name, QWidget* parent = nullptr);

public:
    DECLARE_PLOTTER_GETTER(ImageViewer)

    ImageViewer& setNormalizedModeImage(bool normalize);
    ImageViewer& setAbsoluteModeImage(bool absolute);
    ImageViewer& setFilteredValuesImage(bool filtered);

    ImageViewer& updateUI();
    void displayInfoUnderMouse(const Vector3 &relativeMousePos);

public Q_SLOTS:
    ImageViewer& updateViewOptionsInterface();
};

#endif // IMAGEVIEWER_H
