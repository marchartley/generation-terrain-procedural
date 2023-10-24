#include "Image.h"


Image::Image() : isColor(true), bitDepth(8) {}

Image::Image(const GridF &image)
    : Image()
{
    this->setImage(image);
}

Image::Image(const GridV3 &image)
    : Image()
{
    this->setImage(image);
}


void Image::readColor(png_structp png, png_infop info, png_bytep* row_pointers) {
    int width = png_get_image_width(png, info);
    int height = png_get_image_height(png, info);
    colorImage = GridV3(width, height);

    for (int y = 0; y < height; y++) {
        png_bytep row = row_pointers[y];
        for (int x = 0; x < width; x++) {
            png_bytep px = &(row[x * 4]); // 4 because of RGBA
            colorImage(x, y) = Vector3(px[0], px[1], px[2]);
        }
    }
}

void Image::readBlackWhite(png_structp png, png_infop info, png_bytep* row_pointers) {
    int width = png_get_image_width(png, info);
    int height = png_get_image_height(png, info);
    bwImage = GridF(width, height);

    for (int y = 0; y < height; y++) {
        png_bytep row = row_pointers[y];
        for (int x = 0; x < width; x++) {
            png_bytep px = &(row[x * 4]); // 4 because of RGBA
            bwImage(x, y) = (px[0] + px[1] + px[2]) / 3.0; // A simple average to get grayscale value
        }
    }
}

void Image::writeColor(png_structp png, png_infop info, png_bytep* row_pointers) {
    int width = colorImage.sizeX;
    int height = colorImage.sizeY;

    for (int y = 0; y < height; y++) {
        png_bytep row = row_pointers[y];
        for (int x = 0; x < width; x++) {
            Vector3 px = colorImage(x, y);

            if (bitDepth == 1) {
                uint8_t byteValue = 0;
                if (px.x >= 0.5f) byteValue |= 0b10000000; // R
                if (px.y >= 0.5f) byteValue |= 0b01000000; // G
                if (px.z >= 0.5f) byteValue |= 0b00100000; // B
                row[x / 8] = byteValue; // Here, 8 pixels are packed into a single byte

            } else if (bitDepth == 2) {
                uint8_t r = static_cast<uint8_t>(px.x * 3); // scale 0-1 to 0-3
                uint8_t g = static_cast<uint8_t>(px.y * 3);
                uint8_t b = static_cast<uint8_t>(px.z * 3);
                row[x] = (r << 6) | (g << 4) | (b << 2);

            } else if (bitDepth == 4) {
                uint8_t r = static_cast<uint8_t>(px.x * 15); // scale 0-1 to 0-15
                uint8_t g = static_cast<uint8_t>(px.y * 15);
                uint8_t b = static_cast<uint8_t>(px.z * 15);
                row[x * 2] = (r << 4) | g;
                row[x * 2 + 1] = (b << 4);

            } else if (bitDepth == 8) {
                row[x * 4 + 0] = static_cast<png_byte>(px.x * 255);
                row[x * 4 + 1] = static_cast<png_byte>(px.y * 255);
                row[x * 4 + 2] = static_cast<png_byte>(px.z * 255);
                row[x * 4 + 3] = 255; // Alpha

            } else if (bitDepth == 16) {
                uint16_t r = static_cast<uint16_t>(px.x * 65535);
                uint16_t g = static_cast<uint16_t>(px.y * 65535);
                uint16_t b = static_cast<uint16_t>(px.z * 65535);

                // R channel
                row[x * 8 + 0] = (r >> 8) & 0xFF;
                row[x * 8 + 1] = r & 0xFF;

                // G channel
                row[x * 8 + 2] = (g >> 8) & 0xFF;
                row[x * 8 + 3] = g & 0xFF;

                // B channel
                row[x * 8 + 4] = (b >> 8) & 0xFF;
                row[x * 8 + 5] = b & 0xFF;

                // Alpha channel
                row[x * 8 + 6] = 255;
                row[x * 8 + 7] = 255;
            }
        }
    }
}



void Image::writeBlackWhite(png_structp png, png_infop info, png_bytep* row_pointers) {
    int width = bwImage.sizeX;
    int height = bwImage.sizeY;

    for (int y = 0; y < height; y++) {
        png_bytep row = row_pointers[y];


        if (bitDepth == 1) {
            for (int x = 0; x < width; x += 8) {
                uint8_t byteValue = 0;
                for (int bit = 0; bit < 8 && x + bit < width; ++bit) {
                    float value = bwImage(x + bit, y);
                    if (value >= 0.5f) { // Assuming 0.5 threshold for binary images
                        byteValue |= (1 << (7 - bit));
                    }
                }
                row[x / 8] = byteValue;
            }
        }
        else if (bitDepth == 2) {
            for (int x = 0; x < width; x += 4) {
                uint8_t byteValue = 0;
                for (int bit = 0; bit < 4 && x + bit < width; ++bit) {
                    float value = bwImage(x + bit, y);
                    uint8_t pixelValue = static_cast<uint8_t>(value * 3); // scale 0-1 to 0-3
                    byteValue |= (pixelValue << ((3 - bit) * 2));
                }
                row[x / 4] = byteValue;
            }
        }
        else if (bitDepth == 4) {
            for (int x = 0; x < width; x += 2) {
                float value1 = bwImage(x, y);
                uint8_t pixelValue1 = static_cast<uint8_t>(value1 * 15); // scale 0-1 to 0-15
                uint8_t byteValue = (pixelValue1 << 4);

                if (x + 1 < width) {
                    float value2 = bwImage(x + 1, y);
                    uint8_t pixelValue2 = static_cast<uint8_t>(value2 * 15); // scale 0-1 to 0-15
                    byteValue |= pixelValue2;
                }

                row[x / 2] = byteValue;
            }
        } else {
            for (int x = 0; x < width; x++) {
                float value = bwImage(x, y);
                if (bitDepth == 8) {
                    uint8_t scaledValue = static_cast<uint8_t>(std::round(value * 255)); // Assuming value is in range 0-1
                    row[x] = scaledValue;
                } else if (bitDepth == 16) {
                    uint16_t value16 = static_cast<uint16_t>(std::round(value * 65535)); // Assuming value is in range 0-1
                    row[x * 2] = (value16 >> 8) & 0xFF; // High byte
                    row[x * 2 + 1] = value16 & 0xFF; // Low byte
                }
            }
        }
    }
}

void Image::writeOtherThanPNG(std::string filename, std::string ext)
{
    GridV3 img;
    if (this->isColor) {
        img = this->colorImage;
    } else {
        img = GridV3(this->bwImage.getDimensions());
        for (size_t i = 0; i < img.size(); i++)
            img[i] = Vector3(bwImage[i], bwImage[i], bwImage[i]);
    }
    if (ext.empty())
        ext = getExtension(filename);
    ext = toUpper(ext);

    int nbComp = 4;
    int width = img.sizeX;
    int height = img.sizeY;

    if (ext == "PNG") {
        std::cerr << "Oh oh, I shouldn't be here..." << std::endl;
    } else if (ext == "HDR") {
        std::vector<float> toFloatData(width*height*nbComp);
        for (size_t i = 0; i < img.size(); i++) {
            toFloatData[nbComp * i + 0] = std::max(img[i].x, 0.f);
            toFloatData[nbComp * i + 1] = std::max(img[i].y, 0.f);
            toFloatData[nbComp * i + 2] = std::max(img[i].z, 0.f);
            toFloatData[nbComp * i + 3] = 1.f;
        }
        stbi_write_hdr(filename.c_str(), width, height, nbComp, toFloatData.data());
    } else {

        std::vector<uint8_t> toIntData(width*height*nbComp);

        for (size_t i = 0; i < img.size(); i++) {
            img[i] = Vector3::max(img[i] * 255, Vector3(0, 0, 0));
            toIntData[nbComp * i + 0] = int(img[i].x);
            toIntData[nbComp * i + 1] = int(img[i].y);
            toIntData[nbComp * i + 2] = int(img[i].z);
            toIntData[nbComp * i + 3] = 255;
        }

        if (ext == "JPG")
            stbi_write_jpg(filename.c_str(), width, height, nbComp, toIntData.data(), 95);
        else if (ext == "BMP")
            stbi_write_bmp(filename.c_str(), width, height, nbComp, toIntData.data());
        else if (ext == "TGA")
            stbi_write_tga(filename.c_str(), width, height, nbComp, toIntData.data());
        else {
            throw std::runtime_error("Trying to save image (" + filename + ") without valid extension. Possible extensions :\n\t- png\n\t- jpg\n\t- tga\n\t- bmp\n\t- hdr");
        }
    }
}


Image Image::readFromFile(std::string filename) {
    Image img;

    FILE *fp = fopen(filename.c_str(), "rb");
    if (!fp) {
        throw std::runtime_error("Failed to open the file for reading");
    }

    png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) {
        throw std::runtime_error("Failed to create PNG read struct");
    }

    png_infop info = png_create_info_struct(png);
    if (!info) {
        png_destroy_read_struct(&png, NULL, NULL);
        throw std::runtime_error("Failed to create PNG info struct");
    }

    if (setjmp(png_jmpbuf(png))) {
        png_destroy_read_struct(&png, &info, NULL);
        throw std::runtime_error("Error during PNG read initialization");
    }

    png_init_io(png, fp);
    png_read_info(png, info);

    int width = png_get_image_width(png, info);
    int height = png_get_image_height(png, info);
    png_byte color_type = png_get_color_type(png, info);
    img.bitDepth = png_get_bit_depth(png, info); // Store bitDepth for future use

    if(img.bitDepth == 16)
        png_set_strip_16(png);

    if(color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_palette_to_rgb(png);

    // PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth.
    if(color_type == PNG_COLOR_TYPE_GRAY && img.bitDepth < 8)
        png_set_expand_gray_1_2_4_to_8(png);

    if(png_get_valid(png, info, PNG_INFO_tRNS))
        png_set_tRNS_to_alpha(png);

    // These color_type don't have an alpha channel then fill it with 0xff.
    if(color_type == PNG_COLOR_TYPE_RGB ||
        color_type == PNG_COLOR_TYPE_GRAY ||
        color_type == PNG_COLOR_TYPE_PALETTE)
            png_set_filler(png, 0xFF, PNG_FILLER_AFTER);

    if(color_type == PNG_COLOR_TYPE_GRAY ||
        color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
            png_set_gray_to_rgb(png);

    png_read_update_info(png, info);

    png_bytep* row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * height);
    for (int y = 0; y < height; y++) {
        row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png, info));
    }

    png_read_image(png, row_pointers);

    if (color_type == PNG_COLOR_TYPE_RGB || color_type == PNG_COLOR_TYPE_RGBA) {
        img.readColor(png, info, row_pointers);
    } else if (color_type == PNG_COLOR_TYPE_GRAY) {
        img.readBlackWhite(png, info, row_pointers);
    } else {
        throw std::runtime_error("Unsupported PNG color type");
    }

    for (int y = 0; y < height; y++) {
        free(row_pointers[y]);
    }
    free(row_pointers);

    png_destroy_read_struct(&png, &info, NULL);

    fclose(fp);
    return img;
}

void Image::writeToFile(std::string filename, int desiredBitDepth) {
    std::string extension = toUpper(getExtension(filename));
    if (extension != "PNG")
        return this->writeOtherThanPNG(filename, extension);

    this->bitDepth = desiredBitDepth;
    FILE *fp = fopen(filename.c_str(), "wb");
    if (!fp) {
        throw std::runtime_error("Failed to open the file for writing");
    }

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) {
        fclose(fp);
        throw std::runtime_error("Failed to create PNG write struct");
    }

    png_infop info = png_create_info_struct(png);
    if (!info) {
        png_destroy_write_struct(&png, NULL);
        fclose(fp);
        throw std::runtime_error("Failed to create PNG info struct");
    }

    if (setjmp(png_jmpbuf(png))) {
        png_destroy_write_struct(&png, &info);
        fclose(fp);
        throw std::runtime_error("Error during PNG write initialization");
    }

    png_init_io(png, fp);

    int width, height;
    if (!colorImage.empty()) {
        width = colorImage.sizeX;
        height = colorImage.sizeY;
        png_set_IHDR(png, info, width, height, desiredBitDepth, PNG_COLOR_TYPE_RGBA,
                     PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    } else if (!bwImage.empty()) {
        width = bwImage.sizeX;
        height = bwImage.sizeY;
        png_set_IHDR(png, info, width, height, desiredBitDepth, PNG_COLOR_TYPE_GRAY,
                     PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    } else {
        // Error: Neither colorImage nor bwImage has data
        png_destroy_write_struct(&png, &info);
        fclose(fp);
        throw std::runtime_error("Neither colorImage nor bwImage has data");
    }

    // Write the header information to the file
    png_write_info(png, info);

    png_bytep* row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * height);
    if (!row_pointers) {
        png_destroy_write_struct(&png, &info);
        fclose(fp);
        throw std::runtime_error("Failed to allocate memory for row pointers");
    }

    for (int y = 0; y < height; y++) {
        row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png, info));
        if (!row_pointers[y]) {
            for (int j = 0; j < y; j++) {
                free(row_pointers[j]);
            }
            free(row_pointers);
            png_destroy_write_struct(&png, &info);
            fclose(fp);
            throw std::runtime_error("Failed to allocate memory for a row");
        }
    }

    if (!colorImage.empty()) {
        writeColor(png, info, row_pointers);
    } else if (!bwImage.empty()) {
        writeBlackWhite(png, info, row_pointers);
    }

    png_write_image(png, row_pointers);
    png_write_end(png, NULL);

    for (int y = 0; y < height; y++) {
        free(row_pointers[y]);
    }
    free(row_pointers);

    png_destroy_write_struct(&png, &info);

    fclose(fp);
}

Image &Image::setImage(const GridF &img)
{
    this->colorImage.clear();
    this->bwImage = img;
    this->isColor = false;
    return *this;
}

Image &Image::setImage(const GridV3& img)
{
    this->bwImage.clear();
    this->colorImage = img;
    this->isColor = true;
    return *this;
}


