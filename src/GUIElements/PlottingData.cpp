#include "PlottingData.h"



SingleLineData::SingleLineData(const std::vector<Vector3> &positions, const std::string &name, const Vector3 &color)
    : name(name), color(color)
{
    for (size_t i = 0; i < positions.size(); i++) {
        push_back(SingleLinePoint{positions[i]});
    }
}

void PlotLineData::add(const std::vector<Vector3>& data, const std::string& name, const Vector3& color)
{
    this->push_back(SingleLineData{data, name, color});
    // this->names.push_back(name);
    // this->colors.push_back(color);
}

void PlotLineData::reset()
{
    this->clear();
    // this->names.clear();
    // this->colors.clear();
}




SingleScatterSeries::SingleScatterSeries(const std::vector<Vector3> &positions, const std::vector<std::string> &labels, const Vector3 &color, std::string name)
    : color(color), name(name)
{
    for (size_t i = 0; i < positions.size(); i++) {
        push_back(SingleScatterPoint{positions[i], (labels.size() >= i + 1 ? labels[i] : "")});
    }
}

void PlotScatterData::add(const std::vector<Vector3>& data, const std::vector<std::string>& labels, const Vector3& color, const std::string& series_name)
{
    this->push_back(SingleScatterSeries(data, labels, color, series_name));
    // this->names.push_back(series_name);
    // this->labels.push_back(labels);
    // this->colors.push_back(color);
}

void PlotScatterData::reset()
{
    this->clear();
    // this->names.clear();
    // this->labels.clear();
    // this->colors.clear();
}
