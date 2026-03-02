#ifndef SNAKESEGMENTATIONEDITOR_H
#define SNAKESEGMENTATIONEDITOR_H

#include "GUIElements/ImageViewer.h"
#include "GUIElements/PainterToolsUI.h"

#include "EnvObject/SnakeSegmentation.h"

#include "EnvObject/EnvObject.h"

struct SnakeEditionParameters {
    SnakeSegmentation snake;
    SnakeSegmentationParameters*& params = snake.params;

    /*
    float minConnectivityCost = 0.0f;
    float maxCconnectivityCost  = 1.f;

    float minCurvatureCost = 0.0f;
    float maxCurvatureCost  = 1.f;

    float minImageCost = 0.0f;
    float maxImageCost  = 1.f;

    float minAreaCost = 0.f;
    float maxAreaCost = 1.f;

    float minLengthCost = 0.0f;
    float maxLengthCost  = 1.f;

    float minSlopeCost = 0.f;
    float maxSlopeCost = 1.f;


    float minPositionCost = 0.f;
    float maxPositionCost = 0.f;

    int nbCatapillars = 0;

    float imageBordersCoef = 1.f;
    float imageInsideCoef = 0.f;

    float targetLength = 0.f;
    float targetArea = 0.f;

    bool collapseFirstAndLastPoint = false;*/

};

class SnakeSegmentationEditor : public ImageViewer
{
    Q_OBJECT
protected: // Singleton
    SnakeSegmentationEditor(const std::string& name, QWidget* parent = nullptr);
    SnakeSegmentationEditor(const std::string& name, ChartView* chartView, QWidget* parent = nullptr);

public:
    DECLARE_PLOTTER_GETTER(SnakeSegmentationEditor)

    virtual SnakeSegmentationEditor* updateToolsInterface();
    virtual SnakeSegmentationEditor* updateViewOptionsInterface();

    void showSnakePath();

    SnakeSegmentationEditor* associateEnvObject(EnvObject* obj);

    SnakeEditionParameters snakeParameters;
    PainterToolParams painterParameters;


    GridF currentField;
    GridV3 currentGradientField;

    EnvObject* associatedObject = nullptr;

Q_SIGNALS:
};
#endif // SNAKESEGMENTATIONEDITOR_H
