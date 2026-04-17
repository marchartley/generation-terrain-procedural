#ifndef ENVMATERIALVIEWER_H
#define ENVMATERIALVIEWER_H

#include "GUIElements/ImageViewer.h"
#include "GUIElements/PainterToolsUI.h"

class EnvMaterialViewer : public ImageViewer
{
    Q_OBJECT
public: // protected: // Singleton
    EnvMaterialViewer(const std::string& name, QWidget* parent = nullptr);

public:
    DECLARE_PLOTTER_GETTER(EnvMaterialViewer)

    virtual EnvMaterialViewer& updateToolsInterface();

    PainterToolParams painterParams;

    DECLARE_EVENT(OnImagePainted, (const GridF& newImage), (newImage))
};

#endif // ENVMATERIALVIEWER_H
