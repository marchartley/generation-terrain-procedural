#ifndef FOCUSAREAVIEWER_H
#define FOCUSAREAVIEWER_H

#include "GUIElements/ImageViewer.h"
#include "GUIElements/PainterToolsUI.h"

class FocusAreaViewer : public ImageViewer
{
    Q_OBJECT
protected: // Singleton
    FocusAreaViewer(const std::string& name, QWidget* parent = nullptr);
    FocusAreaViewer(const std::string& name, ChartView* chartView, QWidget* parent = nullptr);

public:
    static FocusAreaViewer* getInstance(std::string name = "");
    static FocusAreaViewer* get(std::string name = "") { return FocusAreaViewer::getInstance(toLower(name)); }
    static FocusAreaViewer* init(const std::string& name, ChartView* chartView = nullptr, QWidget* parent = nullptr);

    virtual FocusAreaViewer* updateToolsInterface();
    virtual FocusAreaViewer* updateViewOptionsInterface();

    // void paintImage(const Vector3 &pos);


    PainterToolParams painterParams;

Q_SIGNALS:
    void imagePainted(const GridF& newImage);
};


#endif // FOCUSAREAVIEWER_H
