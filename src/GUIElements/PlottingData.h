#ifndef PLOTTINGDATA_H
#define PLOTTINGDATA_H

#include <QColor>
#include "GUIElements/PlottingUtils.h"

struct PlotLineData {

    void add(const std::vector<Vector3>& data, const std::string& name, const QColor& color);
    void reset();

    std::vector<std::vector<Vector3>> plot_data;
    std::vector<std::string> plot_names;
    std::vector<QColor> plot_colors;

    bool displayed = true;
};

struct PlotScatterData {

    void add(const std::vector<Vector3>& data, const std::vector<std::string>& labels, const std::vector<QColor>& colors, const std::string& series_name);
    void reset();

    std::vector<std::vector<Vector3>> scatter_data;
    std::vector<std::vector<std::string>> scatter_labels;
    std::vector<std::vector<QColor>> scatter_colors;
    std::vector<std::string> scatter_names;

    std::vector<std::vector<QGraphicsTextItem*>> graphicLabels;

    bool displayed = true;
};

#endif // PLOTTINGDATA_H
