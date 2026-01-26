#include "SnakeSegmentationImplicit.h"

SnakeSegmentationImplicit::SnakeSegmentationImplicit()
    : SnakeSegmentation()
{

}

float SnakeSegmentationImplicit::getImageAt(const Vector3 &p) const
{
    return this->imageField(p);
}

Vector3 SnakeSegmentationImplicit::getGradientImageAt(const Vector3 &p) const
{
    return this->gradientField(p);
}
