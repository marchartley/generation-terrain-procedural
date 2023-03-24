#ifndef PBMREADER_H
#define PBMREADER_H

#include "DataStructure/Matrix3.h"

namespace PbmReader {
Matrix3<float> readGrayscale(std::string path);
Matrix3<Vector3> readColor(std::string path);

int getMagicNumberAndSizes(std::ifstream &file, int& imgW, int& imgH, int &maxValue);

Matrix3<float> readP1(std::ifstream& file, int nW, int nH);
Matrix3<float> readP2(std::ifstream& file, int nW, int nH);
Matrix3<Vector3> readP3(std::ifstream& file, int nW, int nH);
Matrix3<float> readP4(std::ifstream& file, int nW, int nH);
Matrix3<float> readP5(std::ifstream& file, int nW, int nH);
Matrix3<Vector3> readP6(std::ifstream& file, int nW, int nH);
}

#endif // PBMREADER_H
