#ifndef VECTORFIELDDATA_H
#define VECTORFIELDDATA_H

#include "DataStructure/Vector3.h"
#include "DataStructure/Matrix3.h"
#include "DataStructure/Image.h"


struct DisplayedVectorFieldParameters {
    enum VECTOR_DISPLAY { ARROWS, FLOWLINES, NONE };

    VECTOR_DISPLAY displayMode = ARROWS;
    Vector3 backgroundColor = Vector3::white;
    BSpline colorRamp = BSpline({Vector3(70.f, 0.f, 100.f) / 255.f, Vector3(30.f, 160.f, 130.f) / 255.f, Vector3(255.f, 250.f, 0.f)/255.f});
};

struct PlotVectorData {
    PlotVectorData();
    PlotVectorData(const GridV3& field);

    PlotVectorData& setField(const GridV3 &field);

    const GridV3& getField() const { return this->field; }

    std::pair<GridV3, GridF> getFieldImageAndAlpha(const Vector3i &imgSize, const Vector3i &numberOfCells) const;
    static std::pair<GridV3, GridF> createFieldImageAndAlpha(const GridV3& field, Vector3i imgSize, const Vector3i &numberOfCells, DisplayedVectorFieldParameters displayParameters = DisplayedVectorFieldParameters());

    void setOnFieldModified(const std::function<void(void)>& callback) { this->onFieldModifiedCallbacks.push_back(callback); }

    GridV3 field;
    DisplayedVectorFieldParameters displayParameters;

protected:
    std::vector<std::function<void(void)>> onFieldModifiedCallbacks;
    void callOnFieldModifiedCallbacks() { for (const auto& func : onFieldModifiedCallbacks) func(); }
};


#endif // VECTORFIELDDATA_H
