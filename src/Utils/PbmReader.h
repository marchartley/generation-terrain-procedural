#ifndef PBMREADER_H
#define PBMREADER_H

#include "DataStructure/Matrix3.h"

namespace PbmReader {
GridF readGrayscale(std::string path);
GridV3 readColor(std::string path);

int getMagicNumberAndSizes(std::ifstream &file, int& imgW, int& imgH, int &maxValue);

GridF readP1(std::ifstream& file, int nW, int nH);
GridF readP2(std::ifstream& file, int nW, int nH);
GridV3 readP3(std::ifstream& file, int nW, int nH);
GridF readP4(std::ifstream& file, int nW, int nH);
GridF readP5(std::ifstream& file, int nW, int nH);
GridV3 readP6(std::ifstream& file, int nW, int nH);
}

#endif // PBMREADER_H
