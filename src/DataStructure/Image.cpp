#include "Image.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <QSaveFile>
#include <memory>
#include <stdexcept>

Image::Image() : _isColor(true), bitDepth(8) {}

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








void Image::writeOtherThanPNG(const std::string &filename, std::string ext) {
    if (ext.empty())
        ext = toUpper(getExtension(filename));
    const auto dims = _isColor ? colorImage.getDimensions() : bwImage.getDimensions();
    if (dims.x() <= 0 || dims.y() <= 0 || dims.z() != 1)
        throw std::invalid_argument("Image output requires a nonempty 2D grid");
    const int width = dims.x(), height = dims.y();
    const size_t size = size_t(width) * height;
    if (size > size_t(INT_MAX - 128) / 4)
        throw std::length_error("Image exceeds codec byte limits");
    if ((ext == "JPG" || ext == "JPEG" || ext == "TGA") && (width > 65535 || height > 65535))
        throw std::length_error("Image exceeds format dimensions");
    auto sample = [&](size_t i) {
        return _isColor ? colorImage[i] : Vector3(bwImage[i], bwImage[i], bwImage[i]);
    };
    QSaveFile output(QString::fromStdString(filename));
    if (!output.open(QIODevice::WriteOnly))
        throw std::runtime_error("Cannot open image output");
    struct Sink {
        QSaveFile *file;
        bool failed = false;
    } sink{&output};
    auto write = [](void *context, void *bytes, int count) {
        auto *sink = static_cast<Sink *>(context);
        if (sink->file->write(static_cast<const char *>(bytes), count) != count)
            sink->failed = true;
    };
    int ok = 0;
    if (ext == "HDR") {
        std::vector<float> pixels(size * 3);
        for (size_t i = 0; i < size; ++i) {
            auto v = sample(i);
            if (!v.isValid())
                throw std::invalid_argument("Nonfinite HDR sample");
            for (int c = 0; c < 3; ++c)
                pixels[i * 3 + c] = std::max(0.f, v[c]);
        }
        ok = stbi_write_hdr_to_func(write, &sink, width, height, 3, pixels.data());
    } else {
        std::vector<unsigned char> pixels(size * 3);
        for (size_t i = 0; i < size; ++i) {
            auto v = sample(i);
            if (!v.isValid())
                throw std::invalid_argument("Nonfinite image sample");
            for (int c = 0; c < 3; ++c)
                pixels[i * 3 + c] = std::lround(std::clamp(v[c], 0.f, 1.f) * 255);
        }
        if (ext == "JPG" || ext == "JPEG")
            ok = stbi_write_jpg_to_func(write, &sink, width, height, 3, pixels.data(), 95);
        else if (ext == "BMP")
            ok = stbi_write_bmp_to_func(write, &sink, width, height, 3, pixels.data());
        else if (ext == "TGA")
            ok = stbi_write_tga_to_func(write, &sink, width, height, 3, pixels.data());
        else
            throw std::invalid_argument("Unsupported image extension: " + ext);
    }
    if (!ok || sink.failed || !output.commit())
        throw std::runtime_error("Atomic image encoding failed: " + filename);
}


Image Image::readFromFile(const std::string &filename) {
    int width = 0, height = 0, channels = 0;
    if (!stbi_info(filename.c_str(), &width, &height, &channels))
        throw std::runtime_error("Cannot decode image: " + filename);
    Image result;
    result._isColor = channels > 2;
    const int expectedWidth = width, expectedHeight = height;
    auto prepare = [&] {
        if (width != expectedWidth || height != expectedHeight)
            throw std::runtime_error("Image dimensions changed during decoding");
        // A header alone must not force allocation of the much larger float grid.
        if (result._isColor)
            result.colorImage = GridV3(width, height);
        else
            result.bwImage = GridF(width, height);
    };
    auto put = [&](size_t i, float r, float g, float b) {
        if (result._isColor)
            result.colorImage[i] = Vector3(r, g, b);
        else
            result.bwImage[i] = r;
    };
    if (stbi_is_hdr(filename.c_str())) {
        std::unique_ptr<float, decltype(&stbi_image_free)> pixels(
            stbi_loadf(filename.c_str(), &width, &height, &channels, 3), stbi_image_free);
        if (!pixels)
            throw std::runtime_error("Cannot decode HDR image");
        prepare();
        for (size_t i = 0; i < size_t(width) * height; ++i)
            put(i, pixels.get()[3 * i], pixels.get()[3 * i + 1], pixels.get()[3 * i + 2]);
        result.bitDepth = 32;
    } else if (stbi_is_16_bit(filename.c_str())) {
        std::unique_ptr<stbi_us, decltype(&stbi_image_free)> pixels(
            stbi_load_16(filename.c_str(), &width, &height, &channels, 3), stbi_image_free);
        if (!pixels)
            throw std::runtime_error("Cannot decode 16-bit image");
        prepare();
        for (size_t i = 0; i < size_t(width) * height; ++i)
            put(i, pixels.get()[3 * i] / 65535.f, pixels.get()[3 * i + 1] / 65535.f,
                pixels.get()[3 * i + 2] / 65535.f);
        result.bitDepth = 16;
    } else {
        std::unique_ptr<stbi_uc, decltype(&stbi_image_free)> pixels(
            stbi_load(filename.c_str(), &width, &height, &channels, 3), stbi_image_free);
        if (!pixels)
            throw std::runtime_error("Cannot decode image: " + filename);
        prepare();
        for (size_t i = 0; i < size_t(width) * height; ++i)
            put(i, pixels.get()[3 * i] / 255.f, pixels.get()[3 * i + 1] / 255.f,
                pixels.get()[3 * i + 2] / 255.f);
    }
    return result;
}

void Image::writeToFile(const std::string &filename, int depth) {
    const auto dims = _isColor ? colorImage.getDimensions() : bwImage.getDimensions();
    if (dims.x() <= 0 || dims.y() <= 0 || dims.z() != 1)
        throw std::invalid_argument("Image output requires a nonempty 2D grid");
    const int width = dims.x(), height = dims.y(), components = _isColor ? 4 : 1;
    auto value = [&](size_t i, int c) {
        if (_isColor && !colorImage[i].valid)
            throw std::invalid_argument("Invalid RGB pixel");
        float v = _isColor ? (c == 3 ? 1.f : colorImage[i][c]) : bwImage[i];
        if (!std::isfinite(v))
            throw std::invalid_argument("Nonfinite image sample");
        return std::clamp(v, 0.f, 1.f);
    };
    const auto ext = toUpper(getExtension(filename));
    if (ext != "PNG") {
        writeOtherThanPNG(filename, ext);
        return;
    }
    if (depth != 1 && depth != 2 && depth != 4 && depth != 8 && depth != 16)
        throw std::invalid_argument("Unsupported PNG depth");
    if (_isColor && depth < 8)
        throw std::invalid_argument("Color PNG requires 8 or 16 bits");
    const size_t pixelsCount = size_t(width) * height;
    if (pixelsCount > size_t(INT_MAX) / (components * std::max(1, depth / 8)))
        throw std::length_error("PNG staging exceeds supported size");
    const size_t stride = (size_t(width) * components * depth + 7) / 8;
    std::vector<png_byte> pixels(stride * height, 0);
    std::vector<png_bytep> rows(height);
    for (int y = 0; y < height; ++y) {
        rows[y] = pixels.data() + stride * y;
        for (int x = 0; x < width; ++x)
            for (int c = 0; c < components; ++c) {
                auto v = value(size_t(y) * width + x, c);
                size_t sample = size_t(x) * components + c;
                if (depth == 16) {
                    unsigned n = std::lround(v * 65535);
                    rows[y][sample * 2] = n >> 8;
                    rows[y][sample * 2 + 1] = n & 255;
                } else if (depth == 8)
                    rows[y][sample] = std::lround(v * 255);
                else {
                    unsigned n = std::lround(v * ((1 << depth) - 1));
                    rows[y][sample * depth / 8] |= n << (8 - depth - int(sample * depth % 8));
                }
            }
    }
    QSaveFile fp(QString::fromStdString(filename));
    if (!fp.open(QIODevice::WriteOnly))
        throw std::runtime_error("Cannot open image output");
    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    png_infop info = png ? png_create_info_struct(png) : nullptr;
    if (!info) {
        if (png)
            png_destroy_write_struct(&png, nullptr);
        throw std::runtime_error("Cannot create PNG encoder");
    }
    // C resources are released explicitly on libpng's longjmp path. All C++
    // buffers already exist before setjmp; no destructor is skipped by longjmp.
    if (setjmp(png_jmpbuf(png))) {
        png_destroy_write_struct(&png, &info);
        throw std::runtime_error("PNG encoding failed");
    }
    png_set_write_fn(
        png, &fp,
        [](png_structp p, png_bytep bytes, png_size_t count) {
            auto *file = static_cast<QSaveFile *>(png_get_io_ptr(p));
            if (file->write(reinterpret_cast<const char *>(bytes), qint64(count)) != qint64(count))
                png_error(p, "PNG output write failed");
        },
        nullptr);
    png_set_IHDR(png, info, width, height, depth, _isColor ? PNG_COLOR_TYPE_RGBA : PNG_COLOR_TYPE_GRAY,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png, info);
    png_write_image(png, rows.data());
    png_write_end(png, nullptr);
    png_destroy_write_struct(&png, &info);
    if (!fp.commit())
        throw std::runtime_error("Image output close failed");
    bitDepth = depth;
}

const Image &Image::setImage(const GridF &image) {
    colorImage.clear();
    bwImage = image;
    _isColor = false;
    return *this;
}

const Image &Image::setImage(const GridV3 &image) {
    bwImage.clear();
    colorImage = image;
    _isColor = true;
    return *this;
}

GridF Image::getBwImage() const {
    if (!_isColor)
        return bwImage;
    GridF result(colorImage.getDimensions());
    result.iterateParallel(
        [&](size_t i) { result[i] = (colorImage[i].x() + colorImage[i].y() + colorImage[i].z()) / 3.f; });
    return result;
}

GridV3 Image::getColorImage() const {
    if (_isColor)
        return colorImage;
    GridV3 result(bwImage.getDimensions());
    result.iterateParallel([&](size_t i) { result[i] = Vector3(1, 1, 1) * bwImage[i]; });
    return result;
}


