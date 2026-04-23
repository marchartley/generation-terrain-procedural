#ifndef IMAGEDATA_H
#define IMAGEDATA_H

#include "DataStructure/Vector3.h"
#include "DataStructure/Matrix3.h"
#include "DataStructure/Image.h"

#include "Utils/Signals.h"



struct DisplayedImageParameters {
    bool normalized = false;
    bool absolute = false;
    bool clamped = true;
    Vector3 colorRangeMin = Vector3::min();
    Vector3 colorRangeMax = Vector3::max();
    Vector3i displayedColors = Vector3i(1, 1, 1);

    BSpline colorRamp = BSpline({Vector3(1, 0, 0), Vector3(1, 1, 1), Vector3(0, 1, 0)});
};

struct PlotImageData {
    // public:
    PlotImageData();
    PlotImageData(const GridV3& img);
    PlotImageData(const GridF& img);

    PlotImageData& setImage(const GridV3& img);
    PlotImageData& setImage(const GridF& img);
    PlotImageData& setNormalized(bool normalize);
    PlotImageData& setColorRanges(const Vector3& minRange, const Vector3& maxRange);
    PlotImageData& setAbsolute(bool absolute);
    PlotImageData& setClamped(bool clamp);

    GridV3 getImage() const { return this->image.getColorImage(); }
    GridF getImageGrey() const { return this->image.getBwImage(); }
    GridV3 prepareImageForDisplay(const Image &img) const;
    QImage computeDisplayedImage(const Vector3i &imgSize = Vector3i::invalid) const;
    QImage computeDisplayedImage(const GridV3& overlay, const GridF& overlayAlpha) const;
    QImage computeDisplayedImage(const std::map<std::string, GridV3>& overlays,
                                 const std::map<std::string, GridF>& overlayAlphas,
                                 const std::map<std::string, bool>& displayedOverlays,
                                 const std::map<std::string, int>& overlayLayers,
                                 const Vector3i& imgSize) const;

    // void setOnImageModified(const std::function<void(void)>& callback);

    // protected:
    DisplayedImageParameters displayParameters;
    // GridV3 image;
    Image image;

    DECLARE_EVENT(ImageModified, (), ())
// protected:
    // void callOnImageModifiedCallbacks() { for (auto& func : onImageModifiedCallbacks) func(); }
    // std::vector<std::function<void(void)>> onImageModifiedCallbacks;
};


#endif // IMAGEDATA_H
