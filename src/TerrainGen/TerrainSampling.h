#ifndef TERRAINSAMPLING_H
#define TERRAINSAMPLING_H

#include "DataStructure/Matrix3.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

// Scalar fields are sampled at integer coordinates, positive inside, zero on
// the surface. Scale is samples per terrain unit; dimensions are terrain units.
namespace TerrainSampling {
inline Vector3 sampleDimensions(const Vector3& dimensions, const Vector3& scale)
{
    if (!dimensions.isValid() || !scale.isValid() ||
        !std::isfinite(dimensions.x()) || !std::isfinite(dimensions.y()) || !std::isfinite(dimensions.z()) ||
        !std::isfinite(scale.x()) || !std::isfinite(scale.y()) || !std::isfinite(scale.z()) ||
        dimensions.minComp() <= 0.f || scale.minComp() <= 0.f)
        throw std::invalid_argument("Terrain sampling requires positive finite dimensions and scale");
    return (dimensions * scale).roundedUp();
}

inline float heightAt(const GridF& heights, float x, float y)
{
    if (heights.empty()) return 0.f;
    return heights.interpolate(Vector3(std::clamp(x, 0.f, float(heights.sizeX - 1)),
                                        std::clamp(y, 0.f, float(heights.sizeY - 1)), 0));
}

inline GridF heightField(const GridF& heights, const Vector3& dimensions, const Vector3& scale)
{
    GridF result(sampleDimensions(dimensions, scale), -1.f);
    result.iterateParallel([&](int x, int y, int z) {
        result(x, y, z) = heightAt(heights, x / scale.x(), y / scale.y()) - z / scale.z();
    });
    result.defaultValueOnBadCoord = -1.f;
    return result;
}

// Explicit air above the last sample closes even a completely solid column.
inline float aboveColumn(const GridF& values, int x, int y)
{
    return -std::max(std::abs(values(x, y, values.sizeZ - 1)), 1e-6f);
}

inline float topHeight(const GridF& values, int x, int y)
{
    for (int z = int(values.sizeZ) - 1; z >= 0; --z) {
        const float current = values(x, y, z);
        if (current > 0.f) {
            const float above = z + 1 < int(values.sizeZ) ? values(x, y, z + 1) : aboveColumn(values, x, y);
            return z + current / (current - above);
        }
    }
    return 0.f;
}

inline GridF resample(const GridF& source, const Vector3& dimensions, const Vector3& scale)
{
    GridF result(sampleDimensions(dimensions, scale), -1.f);
    auto sampled = source;
    sampled.raiseErrorOnBadCoord = false;
    sampled.defaultValueOnBadCoord = -1.f;
    sampled.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::DEFAULT_VALUE;
    result.iterateParallel([&](int x, int y, int z) {
        result(x, y, z) = sampled.interpolate(Vector3(x, y, z) / scale);
    });
    result.defaultValueOnBadCoord = -1.f;
    return result;
}
}
#endif
