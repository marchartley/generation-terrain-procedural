#ifndef IMAGEPAINTER_H
#define IMAGEPAINTER_H

#include "GUIElements/ImageViewer.h"
#include "GUIElements/PainterToolsUI.h"

class ImagePainter : public ImageViewer
{
    Q_OBJECT
protected: // Singleton
    ImagePainter(const std::string& name, QWidget* parent = nullptr);
    ImagePainter(const std::string& name, ChartView* chartView, QWidget* parent = nullptr);

public:
    static ImagePainter* getInstance(std::string name = "");
    static ImagePainter* get(std::string name = "") { return ImagePainter::getInstance(toLower(name)); }
    static ImagePainter* init(const std::string& name, ChartView* chartView = nullptr, QWidget* parent = nullptr);

    virtual ImagePainter* updateToolsInterface();
    virtual ImagePainter* updateViewOptionsInterface();

    // void paintImage(const Vector3 &pos);


    PainterToolParams painterParams;

Q_SIGNALS:
    void imagePainted(const GridF& newImage);
    void imagePainted(const GridV3& newImage);
};

#endif // IMAGEPAINTER_H
