#ifndef ENVMATERIALVIEWER_H
#define ENVMATERIALVIEWER_H

#include "GUIElements/ImageViewer.h"
#include "GUIElements/PainterToolsUI.h"

class EnvMaterialViewer : public ImageViewer
{
    Q_OBJECT
protected: // Singleton
    EnvMaterialViewer(const std::string& name, QWidget* parent = nullptr);
    EnvMaterialViewer(const std::string& name, ChartView* chartView, QWidget* parent = nullptr);

public:
    // static EnvMaterialViewer* getInstance(std::string name = "");
    // static EnvMaterialViewer* get(std::string name = "") { return EnvMaterialViewer::getInstance(toLower(name)); }
    // static EnvMaterialViewer* init(const std::string& name, ChartView* chartView = nullptr, QWidget* parent = nullptr);
    DECLARE_PLOTTER_GETTER(EnvMaterialViewer)

    virtual EnvMaterialViewer* updateToolsInterface();
    // virtual EnvMaterialViewer* updateViewOptionsInterface();

    // void paintImage(const Vector3 &pos);


    PainterToolParams painterParams;

Q_SIGNALS:
    void imagePainted(const GridF& newImage);
};

#endif // ENVMATERIALVIEWER_H
