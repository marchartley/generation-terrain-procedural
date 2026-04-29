#ifndef PLOTTINGDATA_H
#define PLOTTINGDATA_H

#include "GUIElements/PlottingUtils.h"

#include "Utils/Signals.h"

struct PlotLineData {

    void add(const std::vector<Vector3>& data, const std::string& name, const Vector3& color);
    void reset();

    std::vector<std::vector<Vector3>> data;
    std::vector<std::string> names;
    std::vector<Vector3> colors;

    bool displayed = true;
};

struct PlotScatterData {

    void add(const std::vector<Vector3>& data, const std::vector<std::string>& labels, const Vector3& color, const std::string& series_name);
    void reset();

    std::vector<std::vector<Vector3>> data;
    std::vector<std::vector<std::string>> labels;
    std::vector<Vector3> colors;
    std::vector<std::string> names;

    std::vector<std::vector<QGraphicsTextItem*>> graphicLabels;

    bool displayed = true;
};

#endif // PLOTTINGDATA_H
