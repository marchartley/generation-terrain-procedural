#ifndef IMAGEPAINTER_H
#define IMAGEPAINTER_H

#include "GUIElements/ImageViewer.h"
#include "GUIElements/PainterToolsUI.h"

class ImagePainter : public ImageViewer
{
    Q_OBJECT
public: // protected: // Singleton
    ImagePainter(const std::string& name, QWidget* parent = nullptr);

public:
    // static ImagePainter* getInstance(const std::string& name = "");
    // static ImagePainter* get(const std::string& name = "") { return ImagePainter::getInstance(toLower(name)); }
    // static ImagePainter* init(const std::string& name, std::shared_ptr<ChartView> chartView = nullptr, QWidget* parent = nullptr);
    DECLARE_PLOTTER_GETTER(ImagePainter)

    virtual ImagePainter& updateToolsInterface();
    virtual ImagePainter& updateViewOptionsInterface();

    // void paintImage(const Vector3 &pos);


    PainterToolParams painterParams;

Q_SIGNALS:
    void imagePainted(const GridF& newImage);
    void imagePainted(const GridV3& newImage);
};

#endif // IMAGEPAINTER_H
