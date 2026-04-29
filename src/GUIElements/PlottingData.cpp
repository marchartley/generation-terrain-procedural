#include "PlottingData.h"


void PlotLineData::add(const std::vector<Vector3>& data, const std::string& name, const Vector3& color)
{
    this->data.push_back(data);
    this->names.push_back(name);
    this->colors.push_back(color);
}

void PlotLineData::reset()
{
    this->data.clear();
    this->names.clear();
    this->colors.clear();
}




void PlotScatterData::add(const std::vector<Vector3>& data, const std::vector<std::string>& labels, const Vector3& color, const std::string& series_name)
{
    this->data.push_back(data);
    this->names.push_back(series_name);
    this->labels.push_back(labels);
    this->colors.push_back(color);
}

void PlotScatterData::reset()
{
    this->data.clear();
    this->names.clear();
    this->labels.clear();
    this->colors.clear();
}
