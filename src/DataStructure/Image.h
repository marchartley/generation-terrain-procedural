#ifndef IMAGE_H
#define IMAGE_H

#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "DataStructure/Matrix3.h"

class Image {
private:
    bool isColor;
    int bitDepth; // 8 or 16

    void readColor(png_structp png, png_infop info, png_bytep* row_pointers);
    void readBlackWhite(png_structp png, png_infop info, png_bytep* row_pointers);
    void writeColor(png_structp png, png_infop info, png_bytep* row_pointers);
    void writeBlackWhite(png_structp png, png_infop info, png_bytep* row_pointers);

    void writeOtherThanPNG(std::string filename, std::string ext = "");

public:
    Image();
    Image(const GridF& image);
    Image(const GridV3& image);
    static Image readFromFile(std::string filename);
    void writeToFile(std::string filename, int bit_depth = 16);

    Image& setImage(const GridF& img);
    Image& setImage(const GridV3 &img);

    GridV3 colorImage;
    GridF bwImage;
};

#endif // IMAGE_H
