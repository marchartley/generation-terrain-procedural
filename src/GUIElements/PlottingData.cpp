#include "PlottingData.h"


void PlotLineData::add(const std::vector<Vector3>& data, const std::string& name, const QColor& color)
{
    this->plot_data.push_back(data);
    this->plot_names.push_back(name);
    this->plot_colors.push_back(color);
}

void PlotLineData::reset()
{
    this->plot_data.clear();
    this->plot_names.clear();
    this->plot_colors.clear();
}




void PlotScatterData::add(const std::vector<Vector3>& data, const std::vector<std::string>& labels, const std::vector<QColor>& colors, const std::string& series_name)
{
    this->scatter_data.push_back(data);
    this->scatter_names.push_back(series_name);
    this->scatter_labels.push_back(labels);
    this->scatter_colors.push_back(colors);
}

void PlotScatterData::reset()
{
    this->scatter_data.clear();
    this->scatter_names.clear();
    this->scatter_labels.clear();
    this->scatter_colors.clear();
}
