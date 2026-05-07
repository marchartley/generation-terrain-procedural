#ifndef PLOTTINGDATA_H
#define PLOTTINGDATA_H

#include "GUIElements/PlottingUtils.h"

#include "Utils/Signals.h"

struct SingleLinePoint {
    Vector3 pos;
};

struct SingleLineData : std::vector<SingleLinePoint> {
    SingleLineData(const std::vector<Vector3>& positions, const std::string& name, const Vector3& color);
    std::string name;
    Vector3 color;
    bool displayed = true;
};

struct PlotLineData : public std::vector<SingleLineData> {

    void add(const std::vector<Vector3>& data, const std::string& name, const Vector3& color);
    void reset();

    // std::vector<SingleLineData> data;

    // std::vector<std::vector<Vector3>> data;
    // std::vector<std::string> names;
    // std::vector<Vector3> colors;
    // bool displayed = true;
};

struct SingleScatterPoint {
    Vector3 pos;
    std::string label;
};

struct SingleScatterSeries : public std::vector<SingleScatterPoint> {
    SingleScatterSeries(const std::vector<Vector3>& positions, const std::vector<std::string>& labels, const Vector3& color, std::string name);
    // std::vector<Vector3> data;
    // std::vector<std::string> labels;
    Vector3 color;
    std::string name;
    // std::vector<QGraphicsTextItem*> graphicLabels;

    bool displayed = true;
};

struct PlotScatterData : public std::vector<SingleScatterSeries> {

    void add(const std::vector<Vector3>& data, const std::vector<std::string>& labels, const Vector3& color, const std::string& series_name);
    void reset();

    // std::vector<SingleScatterData> data;
    /*std::vector<std::vector<Vector3>> data;
    std::vector<std::vector<std::string>> labels;
    std::vector<Vector3> colors;
    std::vector<std::string> names;

    std::vector<std::vector<QGraphicsTextItem*>> graphicLabels;

    bool displayed = true;*/
};

#endif // PLOTTINGDATA_H
