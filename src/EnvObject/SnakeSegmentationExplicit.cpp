#include "SnakeSegmentationExplicit.h"

SnakeSegmentationExplicit::SnakeSegmentationExplicit()
    : SnakeSegmentation()
{}



float SnakeSegmentationExplicit::getImageAt(const Vector3 &p) const
{
    return this->image.interpolate(p);
}

Vector3 SnakeSegmentationExplicit::getGradientImageAt(const Vector3 &p) const
{
    return this->gradientField.interpolate(p);
}
