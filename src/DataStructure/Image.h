#ifndef IMAGE_H
#define IMAGE_H

#include <png.h>
#include <stdlib.h>
#include "DataStructure/Matrix3.h"


#include "Utils/stb_image.h"
#include "Utils/stb_image_write.h"
#include "Utils/stb_image_resize.h"

class Image {
private:
    bool _isColor;
    int bitDepth; // 8 or 16

    void writeOtherThanPNG(const std::string& filename, std::string ext = "");

public:
    Image();
    Image(const GridF& image);
    Image(const GridV3& image);
    static Image readFromFile(const std::string& filename);
    void writeToFile(const std::string& filename, int bit_depth = 16);

    const Image& setImage(const GridF& img);
    const Image &setImage(const GridV3 &img);

    GridF getBwImage() const;
    GridV3 getColorImage() const;

    bool isColor() const { return this->_isColor; }

    bool empty() const { return this->colorImage.empty() && this->bwImage.empty(); }

    GridV3 colorImage;
    GridF bwImage;
};

#endif // IMAGE_H
