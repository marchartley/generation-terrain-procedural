#ifndef SNAKESEGMENTATIONEXPLICIT_H
#define SNAKESEGMENTATIONEXPLICIT_H

#include "EnvObject/SnakeSegmentation.h"

class SnakeSegmentationExplicit : public SnakeSegmentation
{
public:
    SnakeSegmentationExplicit();


    float getImageAt(const Vector3& p) const;
    Vector3 getGradientImageAt(const Vector3& p) const;


    GridF image;         // Grayscale image grid
    GridV3 gradientField; // Gradient field of the image
};

#endif // SNAKESEGMENTATIONEXPLICIT_H
