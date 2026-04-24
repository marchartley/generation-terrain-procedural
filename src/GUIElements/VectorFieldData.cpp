#include "VectorFieldData.h"

#include "GUIElements/PlottingUtils.h"


PlotVectorData::PlotVectorData() : PlotVectorData(GridV3())
{}

PlotVectorData::PlotVectorData(const GridV3 &field) : field(field)
{}

PlotVectorData& PlotVectorData::setField(const GridV3 &field)
{
    this->field = field;
    // this->callOnFieldModifiedCallbacks();
    emitFieldModified();
    return *this;
}

std::pair<GridV3, GridF> PlotVectorData::getFieldImageAndAlpha(const Vector3i &imgSize, const Vector3i& numberOfCells) const
{
    return PlotVectorData::createFieldImageAndAlpha(this->field, imgSize, numberOfCells, this->displayParameters);
}

std::pair<GridV3, GridF> PlotVectorData::createFieldImageAndAlpha(const GridV3 &field, Vector3i imgSize, const Vector3i &numberOfCells, DisplayedVectorFieldParameters displayParameters)
{
    const Vector3& backgroundColor = displayParameters.backgroundColor;
    if (!imgSize.isValid()) imgSize = field.getDimensions();
    GridV3 img(imgSize.x(), imgSize.y(), 1, backgroundColor);
    GridF alpha(img.getDimensions());
    if (field.empty()) return {img, alpha};

    GridV3 reduced = field.resize(numberOfCells);

    float minMag = std::numeric_limits<float>::max();
    float maxMag = std::numeric_limits<float>::lowest();
    reduced.iterate([&] (size_t i) {
        float mag = reduced[i].norm2();
        minMag = std::min(minMag, mag);
        maxMag = std::max(maxMag, mag);
    });
    minMag = std::sqrt(minMag);
    maxMag = std::sqrt(maxMag);
    // std::cout << minMag << " " << maxMag << std::endl;
    Vector3 reducedSize = reduced.getDimensions(); //imgSize / numberOfCells;
    Vector3 imageToReducedRatio = Vector3((float)imgSize.x() / (float)reducedSize.x(), (float)imgSize.y() / (float)reducedSize.y(), 1.f);
    Vector3 fieldToReducedRatio = Vector3((float)field.sizeX / (float)reducedSize.x(), (float)field.sizeY / (float)reducedSize.y(), 1.f);
    Vector3 fieldToImageRatio = Vector3((float)field.sizeX / (float)imgSize.x(), (float)field.sizeY / (float)imgSize.y(), 1.f);

    if (displayParameters.displayMode == DisplayedVectorFieldParameters::ARROWS) {
        reduced.iterateParallel([&] (const Vector3& _p) {
            Vector3 p = _p + Vector3(.5f, .5f);
            // AABBox cell((p - Vector3(.5f, .5f, 1)) * ratio, (p + Vector3(.5f, .5f, 1)) * ratio); // Added an depth (z) to avoid issue on the intersection computation
            Vector3 vec = reduced.interpolate(p);
            if (!vec.isValid()) return;
            float mag = vec.norm();
            if (mag < 1e-5) return;
            Vector3 dir = vec / mag;
            Vector3 color = displayParameters.colorRamp.getPoint(0.f);
            float relativeMag = 1.f;
            if (std::abs(minMag - maxMag) > 1e-5) {
                relativeMag = interpolation::linear(mag, minMag, maxMag);
                color = colorPalette(relativeMag, displayParameters.colorRamp.getPath());
            }
            bool valid = dir.xy().norm2() > 1e-5;
            if (!valid) return;

            Vector3 startLine = (p - dir * interpolation::inv_linear(relativeMag, .5f, 1.f)) * imageToReducedRatio;
            Vector3 endLine = (p + dir * interpolation::inv_linear(relativeMag, .5f, 1.f))  * imageToReducedRatio;
            float length = (endLine - startLine).norm();

            img = PlottingUtils::drawLine(img, color, startLine, endLine);
            alpha = PlottingUtils::drawLine(alpha, 1.f, startLine, endLine);

            img = PlottingUtils::drawLine(img, color, endLine, endLine - dir.rotated(deg2rad(20), Vector3(0, 0, 1)) * length * .3f);
            alpha = PlottingUtils::drawLine(alpha, 1.f, endLine, endLine - dir.rotated(deg2rad(20), Vector3(0, 0, 1)) * length * .3f);

            img = PlottingUtils::drawLine(img, color, endLine, endLine - dir.rotated(deg2rad(-20), Vector3(0, 0, 1)) * length * .3f);
            alpha = PlottingUtils::drawLine(alpha, 1.f, endLine, endLine - dir.rotated(deg2rad(-20), Vector3(0, 0, 1)) * length * .3f);
        });
    }
    else if (displayParameters.displayMode == DisplayedVectorFieldParameters::FLOWLINES) {
        int trailLength = 100;
        float stepLength = .1f;

        for (int x = 0; x < numberOfCells.x() - 1; x++) {
            for (int y = 0; y < numberOfCells.y() - 1; y++) {
                Vector3 p = Vector3(x + .5f, y + .5f) * fieldToReducedRatio;
                for (int i = 0; i < trailLength; i++) {
                    Vector3 dir = field.interpolate(p);
                    if (!dir.isValid()) break;
                    float mag = dir.norm();
                    if (mag < 1e-5) break;
                    auto color = colorPalette(interpolation::inv_linear(mag, minMag, maxMag), displayParameters.colorRamp.getPath());

                    dir = (dir.maxMagnitude(1.f)) * fieldToReducedRatio * stepLength;
                    Vector3 end = p + dir;

                    img = PlottingUtils::drawLine(img, color, p / fieldToImageRatio, end / fieldToImageRatio);
                    alpha = PlottingUtils::drawLine(alpha, 1.f, p / fieldToImageRatio, end / fieldToImageRatio);

                    p = end;
                    if (x == 0 && y == 0) std::cout << i << " -> " << p << "(" << mag << " / " << minMag << " / " << maxMag << ")" << std::endl;
                }
            }
        }
    }
    return {img, alpha};
}
