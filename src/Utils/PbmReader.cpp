#include "PbmReader.h"


std::string readLineAndIgnoreComments(std::ifstream& file) {
    std::string line;
    while(std::getline(file, line)) {
        if (startsWith(line, "#")) {
            continue;
        }
        return line;
    }
    return "";
}

std::vector<float> extractAllValuesFromLine(std::string line, bool isBinary) {
    std::vector<float> values;
    if (isBinary) {
//        std::vector<std::string> items = split(line);
        std::vector<unsigned char> items(line.begin(), line.end());
        for (const auto& item : items) {
            values.push_back((int)item);
        }
    } else {
        std::vector<std::string> items = split(line, " ");
        for (const auto& item : items) {
            values.push_back(std::stof(item));
        }
    }
    return values;
}

std::vector<float> extractAllValuesFromFile(std::ifstream& file, bool isBinary) {
    std::vector<float> values;
    std::string line;
    while (!(line = readLineAndIgnoreComments(file)).empty()) {
        std::vector<float> newValues = extractAllValuesFromLine(line, isBinary);
        values.insert(values.end(), newValues.begin(), newValues.end());
    }
    return values;
}

GridF PbmReader::readGrayscale(std::string path)
{
    std::ifstream file(path);
    int magicNumber, nW, nH, maxValue;
    magicNumber = getMagicNumberAndSizes(file, nW, nH, maxValue);
    if (magicNumber == 1) return readP1(file, nW, nH) / (float)maxValue;
    else if (magicNumber == 2) return readP2(file, nW, nH) / (float)maxValue;
    else if (magicNumber == 4) return readP4(file, nW, nH) / (float)maxValue;
    else if (magicNumber == 5) return readP5(file, nW, nH) / (float)maxValue;
    else {
        std::cerr << "Unable to read grayscale image '" << path << "'. Maybe it's a color image..." << std::endl;
        return {};
    }
}

GridV3 PbmReader::readColor(std::string path)
{
    std::ifstream file(path);
    int magicNumber, nW, nH, maxValue;
    magicNumber = getMagicNumberAndSizes(file, nW, nH, maxValue);
    if (magicNumber == 3) return readP3(file, nW, nH) / maxValue;
    else if (magicNumber == 6) return readP6(file, nW, nH) / maxValue;
    else {
        std::cerr << "Unable to read color image '" << path << "'. Maybe it's a grayscale image..." << std::endl;
        return {};
    }
}

int PbmReader::getMagicNumberAndSizes(std::ifstream& file, int& imgW, int& imgH, int& maxValue)
{
    std::string line;
    imgW = 0;
    imgH = 0;
    int number = -1;
    while (!(line = readLineAndIgnoreComments(file)).empty()) {
        if (startsWith(line, "P")) {
            number = std::stoi(line.substr(line.find("P") + 1));
        }
        else if (number > 0 && imgW == 0 && imgH == 0) {
            auto values = extractAllValuesFromLine(line, false);
            if (values.size() > 0) {
                imgW = values[0];
                imgH = values[1];
            }
        }
        else if (imgW > 0 && imgH > 0) {
            auto values = extractAllValuesFromLine(line, false);
            if (values.size() > 0) {
                maxValue = values[0];
                return number;
            }
        }
    }

    //std::cerr << "Did not find PBM magic number (P[1-6]) in file '" << file << "'. There must be an error somewhere wuth this file." << std::endl;
    return -1;
}

GridF PbmReader::readP1(std::ifstream& file, int nW, int nH)
{
    GridF img(nW, nH);
    auto values = extractAllValuesFromFile(file, false);
    img.data = values;
    return img;
}

GridF PbmReader::readP2(std::ifstream& file, int nW, int nH)
{
    GridF img(nW, nH);
    auto values = extractAllValuesFromFile(file, false);
    img.data = values;
    return img;
}

GridV3 PbmReader::readP3(std::ifstream& file, int nW, int nH)
{
    GridV3 img(nW, nH);
    auto values = extractAllValuesFromFile(file, false);
    for (size_t i = 0; i < values.size() / 3; i++)
        img[i] = Vector3(
                    values[i*3 + 0],
                    values[i*3 + 0],
                    values[i*3 + 0]
                );
    return img;
}

GridF PbmReader::readP4(std::ifstream& file, int nW, int nH)
{
    GridF img(nW, nH);
    auto values = extractAllValuesFromFile(file, true);
    img.data = values;
    return img;
}

GridF PbmReader::readP5(std::ifstream& file, int nW, int nH)
{
    GridF img(nW, nH);
    auto values = extractAllValuesFromFile(file, true);
    img.data = values;
    return img;
}

GridV3 PbmReader::readP6(std::ifstream& file, int nW, int nH)
{
    GridV3 img(nW, nH);
    auto values = extractAllValuesFromFile(file, true);
    for (size_t i = 0; i < values.size() / 3; i++)
        img[i] = Vector3(
                    values[i*3 + 0],
                    values[i*3 + 0],
                    values[i*3 + 0]
                );
    return img;
}
