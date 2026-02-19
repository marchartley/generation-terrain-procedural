#ifndef FANCYSLIDER_H
#define FANCYSLIDER_H

#include <QObject>
#include <QSlider>

class FancySlider : public QSlider{
    Q_OBJECT
public:
    explicit FancySlider(float multiplier = 1.0, QWidget *parent = nullptr);
    explicit FancySlider(Qt::Orientation orientation, float min = 0.0, float max = 100.0, float multiplier = 1.0, QWidget *parent = nullptr);
    virtual ~FancySlider();

    float getfValue();
    void setfValue(float val);
    void setfRange(float min, float max);
    void addTicks(const std::vector<std::pair<float, std::string>>& posAndLabels);
    float multiplier;

    std::vector<std::pair<float, std::string>> ticksPositionsAndLabels;

    void paintEvent(QPaintEvent *ev);

Q_SIGNALS:
    void floatValueChanged(float value);
    void doubleClicked();

public Q_SLOTS:
    void notifyValueChanged(int value);
    void mouseDoubleClickEvent(QMouseEvent *event) { Q_EMIT this->doubleClicked(); }
protected:
    virtual void sliderChange(SliderChange change);
};

#endif // FANCYSLIDER_H
