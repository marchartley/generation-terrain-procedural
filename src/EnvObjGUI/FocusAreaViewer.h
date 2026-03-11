#ifndef FOCUSAREAVIEWER_H
#define FOCUSAREAVIEWER_H

#include "GUIElements/ImageViewer.h"
#include "GUIElements/PainterToolsUI.h"

class FocusAreaViewer : public ImageViewer
{
    Q_OBJECT
public: // protected: // Singleton
    FocusAreaViewer(const std::string& name, QWidget* parent = nullptr);

public:
    DECLARE_PLOTTER_GETTER(FocusAreaViewer)

    virtual FocusAreaViewer& updateToolsInterface();
    PainterToolParams painterParams;

Q_SIGNALS:
    void imagePainted(const GridF& newImage);
};


#endif // FOCUSAREAVIEWER_H
