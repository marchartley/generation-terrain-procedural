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

    void setfValue(float val);
    void setfRange(float min, float max);
    float multiplier;

Q_SIGNALS:
    void floatValueChanged(float value);

public Q_SLOTS:
    void notifyValueChanged(int value);
protected:
    virtual void sliderChange(SliderChange change);
};

#endif // FANCYSLIDER_H
